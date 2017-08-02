/*
 * PHLAWD: Phylogeny assembly with databases
 * Copyright (C) 2010  Stephen A. Smith
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

/*
 * SQLiteConstructor.cpp
 */

#include "SWPS3_matrix.h"

#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <limits>
#include <set>

using namespace std;

#include "omp.h"
#include "libsqlitewrapped.h"
#include "tree.h"
#include "tree_reader.h"
#include "sequence.h"
#include "fasta_util.h"
#include "SQLiteConstructor.h"
#include "SQLiteDBController.h"
#include "utils.h"

//the number for the mafft alignment sample
#define RANDNUM 1000
#define MAX_SEQS_PER_ALIGNMENT 7000		// TODO: make these user-settable
#define SHRINKABLE_THRESHOLD 10000
#define MAD_CONSTANT 1.4826
#define MAX_UPDATE_PROPORTION 0.2

//public

template<class T>
inline std::string to_string(const T& t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

struct seqRecordTuple {
	string db_record_id;
	string ncbi_tax_id;
};

SQLiteConstructor::SQLiteConstructor(
		string cn,
		vector<string> searchstr,
		bool searchlit,
		bool exclsearchlit,
		string genen,
		string genedb,
		double mad_cut,
		double cover,
		double ident,
		string dbs,
		string known_seq_filen,
		bool its,
		int numt,
		bool autom,
		bool inupdatedb,
		string inupdatefile,
		bool assignleft,
		int shrinkthresh,
        bool labeluserseqs):
		clade_name(cn),
				search(searchstr),
				searchliteral(searchlit),
				exclsearchliteral(exclsearchlit),
				gene_name(genen),
				gene_db_name(genedb),
				mad_cutoff(mad_cut),
				coverage(cover),
				identity(ident),
				db(dbs),
				useITS(its),
				numthreads(numt),
				automated(autom),
				updateDB(inupdatedb),
				userfasta(false),
				assignleftovers(assignleft),
				shrinkablethreshold(shrinkthresh),
                labelUserSequences(labeluserseqs){

	FastaUtil seqReader;
	//added updating seqs from db
	//added updating seqs from file
	if (inupdatefile.length() > 0) {
		updateFILE = true;
		updatef = inupdatefile;
	} else {
		updateFILE = false;
	}
	//initialize gene database
	gene_db = GeneDB(gene_db_name);
	
	// open the source db
	SQLiteDBController dbc = SQLiteDBController(db);
	
	genefoldername = genen + "_TEMPFILES/";
	known_seqs = new vector<Sequence>();
	user_seqs = new vector<Sequence>();
	seqReader.read_user_fasta_into(*known_seqs, known_seq_filen); // the keepfile sequences cannot have gaps
	ncbi_saturation = true;
	userskipsearchdb = false;
	skipdbcheck = false;
	justseqquery = false;
}

/*
 * This will supersede the previous one in the main file. If you have higher taxa
 * in the listfile and you say they are wildcards then it will supersede if you
 * have higher taxa and don't want to pick one as above.
 */
void SQLiteConstructor::set_only_names_from_file(string filename, bool containshi, bool containswi) {
	onlynamesfromfile = true;
	containshigher = containshi;
	containswild = containswi;
	if (containswild == true)
		containshigher = false;
	listfilename = filename;
}

void SQLiteConstructor::set_exclude_names_from_file(string filename, bool containshi, bool containswi) {
	excludenamesfromfile = true;
	excludefilename = filename;
	containshigherex = containshi;
	containswildex = containswi;
	if (containswildex == true)
		containshigherex = false;
}

void SQLiteConstructor::set_exclude_gi_from_file(string filename) {
	excludegifromfile = true;
	exclude_gi_filename = filename;
}

void SQLiteConstructor::set_include_gi_from_file(string filename) {
	includegifromfile = true;
	include_gi_filename = filename;
}

void SQLiteConstructor::set_user_skip_search() {
	userskipsearchdb = true;
}

void SQLiteConstructor::set_justseqquery(bool setit) {
	justseqquery = setit;
}

/*
 * if the user guide tree covers less than 50% of the taxa, falling back
 * on NCBI
 * this will allow for skipping the check to the ncbi db, not really that 
 * useful, but good for things like simulated data
 */
void SQLiteConstructor::set_user_guide_tree(string filename, bool inskipdbcheck) {
	usertree = true;
	usertreefile = filename;
	//read the tree here
	TreeReader nw;
	ifstream infile2(filename.c_str());
	if (!infile2) {
		cerr << "user guide tree: " << filename << " doesn't exist" << endl;
		exit(0);
	}
	vector<string> lines;
	string line;
	while (getline(infile2, line)) {
		lines.push_back(line);
	}
	infile2.close();
	userguidetree = nw.readTree(lines[0]);
	cout << "read user guide tree with : " << userguidetree->getExternalNodeCount() << " tips" << endl;
	//need to number the internal nodes so that they can be used for the run
	int count = 0;
	for (int i = 0; i < userguidetree->getInternalNodeCount(); i++) {
		userguidetree->getInternalNode(i)->setName(to_string(count));
		count += 1;
	}
	//procedure would be to
	//get the ncbi taxa for each of these
	//associating that information with each of the tip nodes as comment
	skipdbcheck = inskipdbcheck;
	if (skipdbcheck == false) {
		Database conn(db);
		//TODO: need to make this faster
		cout << "matching user guide tree names: " << userguidetree->getExternalNodeCount() << endl;
		for (int i = 0; i < userguidetree->getExternalNodeCount(); i++) {
			string tname = userguidetree->getExternalNode(i)->getName();
			//search for the name in the database
			string sql = "SELECT ncbi_id FROM taxonomy WHERE edited_name = '" + tname + "' OR ncbi_id = '" + tname + "';";
			Query query(conn);
			query.get_result(sql);
			int count1 = 0;
			bool nameset = false;
			string nameval;
			while (query.fetch_row()) {
				nameset = true;
				nameval = to_string(query.getval());
				count1 += 1;
			}
			query.free_result();
			if (nameset == true) {
				userguidetree->getExternalNode(i)->setComment(nameval);
				//cout << tname << "="<<nameval<<endl;
			} else {
				cerr << tname << " is not in the ncbi database as a number or name" << endl;
			}
		}
	}
	//this will change the search to use the user tree instead of the ncbi tree
	ncbi_saturation = false;
}

/*
 * allows for adding user seqs
 * will also allow for skipping a db check. really not that useful unless
 * there is simulated data or something similar
 */
void SQLiteConstructor::set_user_fasta_file(string filename, bool skipdbcheck) {

	// TODO: need to test if this works with the new sequence objects

	userfasta = true;
	userfastafile = filename;
	FastaUtil fu;
	cout << "reading user fasta file: " << userfastafile << endl;
	fu.read_user_fasta_into(*user_seqs, userfastafile);
	cout << "successfully read " << user_seqs->size() << " user seqs" << endl;

	//need to edit the names if there aren't already changed
	//should go into the log
	/*	for (int i = 0; i < user_seqs->size(); i++) {
	 string tstring(user_seqs->at(i).get_label());
	 fix_bad_chars_for_seq_names(tstring);

	 //if user is not in the front then add it

	 if (tstring.find("user_") != 0)
	 tstring = "user_" + tstring;
	 cout << "changing " << user_seqs->at(i).get_label() << " -> " << tstring << endl;

	 // TODO: not sure what this is getting set as
	 //		user_seqs->at(i).set_label(tstring);
	 } */

	//if the ids can be found in the database, this will go ahead and make that link
	//this should be able to link to the genedb eventually as well
	//store the ncbi id in the comments
	if (skipdbcheck == false) {
		Database conn(db);
		//TODO: need to make this faster
		cout << "trying to match sequence names to ncbi taxon names " << endl;

		set<string> matched_ncbi_ids;
		set<string> observed_unmatched_names;

		for (int i = 0; i < user_seqs->size(); i++) {
//			string tname = user_seqs->at(i).get_label().substr(5, user_seqs->at(i).get_label().size());

			string tname = user_seqs->at(i).get_taxon_name();

			//search for the name in the database
			string sql = "SELECT ncbi_id FROM taxonomy WHERE edited_name = '" + tname + "' OR ncbi_id = '" + tname + "';";
			Query query(conn);
			query.get_result(sql);
			int count1 = 0;
			bool nameset = false;
			string nameval;
			while (query.fetch_row()) {
				nameset = true;
				nameval = to_string(query.getval());
				count1 += 1;
			}
			query.free_result();

			// if we found a taxon, record it
			if (nameset == true) {

				// but don't allow multiple sequences for the same taxon
				if (matched_ncbi_ids.find(nameval) == matched_ncbi_ids.end()) {
					user_seqs->at(i).set_ncbi_tax_id(nameval);
					cout << tname << "=" << nameval << endl;
					matched_ncbi_ids.insert(nameval);
				} else {
					cerr << endl << "No more than one user sequence may be defined for any ncbi taxon, but there are at least two sequences matching: " << tname << " = " << nameval << "." << endl << endl;
					exit(0);
				}

			// didn't find an ncbi taxon
			} else {

				// don't allow identical sequence names
				if (observed_unmatched_names.find(tname) == observed_unmatched_names.end()) {
					cout << tname << " is not in the ncbi database as a number or name" << endl;
					observed_unmatched_names.insert(tname);
				} else {
					cerr << endl << "User sequence names must be unique. There are at least two sequences named: " << tname << "." << endl << endl;
					exit(0);
				}
			}
            
            if (labelUserSequences) {
                // will ensure the name doesn't conflict with other sequences for this taxon
                // but there will be duplicates and user seqs won't be chimera-izable with gb seqs from other alignments
                user_seqs->at(i).set_taxon_name(tname + "_usersequence");
            }
		}
	}
	if (usertree == true) {
		cout << "matching user fasta seqs to user guide tree" << endl;
		int count = 0;
		for (int i = 0; i < userguidetree->getExternalNodeCount(); i++) {
			string tname = userguidetree->getExternalNode(i)->getName();
			// cout << tname << endl;
			for (int j = 0; j < user_seqs->size(); j++) {
				if (tname == user_seqs->at(j).get_taxon_name() || tname == user_seqs->at(j).get_comment()
						|| ("user_" + tname == user_seqs->at(j).get_taxon_name())
						|| ("user_" + tname == user_seqs->at(j).get_comment())) {
					user_fasta_node_map[&(user_seqs->at(j))] = userguidetree->getExternalNode(i);
					count += 1;
				}
			}
		}
		cout << "matches: " << count << " prop:" << count / user_seqs->size() << endl;
	}
	cout << "finished reading user fasta file" << endl;
}

int SQLiteConstructor::run() {

	string logn = gene_name + ".log";
	logfile.open(logn.c_str());
	string gin = gene_name + ".gi";
	string user_fasta_fname = gene_name + ".userfasta";

	write_EDNAFILE();    // if it doesn't exist

	// if we need to make a new gene db, do it
	bool overwrite_gene_db;
	if (updateDB == true)
		overwrite_gene_db = false;
	else
		overwrite_gene_db = true;
	gene_db.initialize(overwrite_gene_db);

	// set up containers in case we're doing an update
	vector<Sequence> preexisting_seqs;
	vector<string> preexisting_db_seqs_ncbi_ids;
	vector<string> preexisting_user_seqs_names;

	if (updateDB == true) {
		cout << "processing existing sequences" << endl;

		// prepare for update by reading in all the preexisting sequences
		gene_db.load_all_original_sequences_into(preexisting_seqs);
		for (int i = 0; i < preexisting_seqs.size(); i++) {

			// sort preexisting sequences by whether they are ncbi seqs or user-supplied
			if (preexisting_seqs[i].is_user_seq())
				preexisting_user_seqs_names.push_back(preexisting_seqs[i].get_taxon_name());
			else
				preexisting_db_seqs_ncbi_ids.push_back(preexisting_seqs[i].get_ncbi_tax_id());
		}
		cout << "processed " << preexisting_seqs.size() << " pre-existing seqs, including " << preexisting_db_seqs_ncbi_ids.size() << " of ncbi origin and "
				<< preexisting_user_seqs_names.size() << " user-supplied ones" << endl;

	} else { // not an update run

		gifile.open(gin.c_str(), fstream::out);
		gifile << "ncbi_tax_id\tgi_number\tedited_name" << endl;

		// output the user fasta seqs
		if (userfasta == true) {
			ufafile.open(user_fasta_fname.c_str(), fstream::out);
			ufafile << "edited_name\tncbi_taxon_id" << endl;
		}
	}

	// create a folder to hold temporary files for this run
	mkdir(genefoldername.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);

	// set up containers to hold seqs we find in the source db
	vector<seqRecordTuple> all_seqs_matching_search;
	vector<Sequence> filtered_starting_seqs;

	// prepare to search source db
	Database conn(db);
	int search_clade_id;
	string s_search_clade_id;

	if (userskipsearchdb == false) {

		// find *all* the sequences in the source db that satisfy the search criteria
		load_info_for_all_seqs_matching_search_into(all_seqs_matching_search);

		if (all_seqs_matching_search.size() == 0) {
			cerr << "there are no seqs in the source db that match the search criteria. quitting";
			exit(0);
		}

		cout << "connected to " << db << endl;
		vector<int> found_ids;

		// get the ncbi taxon id for the search clade
		if (automated == true) {
			// for so-called 'automated' runs, we use an ncbi id instead of a name?
			// this seems redundant... if we already know the ncbi id, why are we attempting to select it? just for validation?
			// TODO: this is going to be broken now.
			Query query(conn);
			query.get_result("SELECT ncbi_id FROM taxonomy WHERE ncbi_id = " + clade_name);
			while (query.fetch_row()) {
				found_ids.push_back(query.getval());
			}
			query.free_result();

		} else {
			// for normal runs, locate the the search taxon by name
			Query query(conn);
			query.get_result("SELECT ncbi_id FROM taxonomy WHERE name == '" + clade_name + "' and name_class == 'scientific name'");
			while (query.fetch_row()) {
				found_ids.push_back(query.getval());
			}
			query.free_result();
		}

		// print all the ncbi taxon ids that match the search clade
		cout << "found " << found_ids.size() << " ncbi taxon ids matching " << clade_name << ":" << endl;
		for (int i = 0; i < found_ids.size(); ++i) {
			cout << found_ids[i] << endl;
		}

		// assume the first matching id is the clade we want
		search_clade_id = found_ids[0];
		s_search_clade_id = to_string(found_ids[0]).c_str();
		cout << "will only be using the first one: " << search_clade_id << endl;

		// filter out any sequences that matched the search criteria but are not within the search clade
		filtered_starting_seqs = make_seqs_from_seq_tuples_for_taxon(search_clade_id, all_seqs_matching_search);
		cout << "\nfound " << filtered_starting_seqs.size() << " matching sequences in the source db" << endl << endl;

		if (excludegifromfile == true) { // if excluding gi's from file, filter against those
			filtered_starting_seqs = exclude_gis_from_file(filtered_starting_seqs);
			cout << "after excluding gis in the gis to exclude file: " << filtered_starting_seqs.size() << endl << endl;
		}

		if (includegifromfile == true) { // if including gi's from file, add those
			filtered_starting_seqs = include_gis_from_file(filtered_starting_seqs);
			cout << "after including gis in the gis to include file: " << filtered_starting_seqs.size() << endl << endl;
		}

		if (onlynamesfromfile == true) { // if only using the names from a list, filter against the list
			filtered_starting_seqs = use_only_names_from_file(filtered_starting_seqs);
			cout << "after filtering taxa not in the include file: " << filtered_starting_seqs.size() << endl << endl;
		}

		if (excludenamesfromfile == true) { // if excluding names from file, filter them out
			filtered_starting_seqs = exclude_names_from_file(filtered_starting_seqs);
			cout << "after filtering against search terms in the exclude file: " << filtered_starting_seqs.size() << endl << endl;
		}

		if (filtered_starting_seqs.size() == 0) {
			cout << "there were no seqs left after the filtering steps. phlawd will exit" << endl;
			exit(0);
		}
	}    // end userskipsearchdb == false

	reduce_genomes(&filtered_starting_seqs);

	// initialize container to hold the seqs that pass the coverage/identity cutoffs
	vector<Sequence> * seqs_to_align = new vector<Sequence>();

	// run the smith-waterman seq similiarity tests
	cout << "running similarity tests (coverage/identity) against keepfile" << endl;
	if (justseqquery == true)
		get_best_hits_openmp_SWPS3_justquery(filtered_starting_seqs, seqs_to_align);
	else
		get_best_hits_openmp_SWPS3(filtered_starting_seqs, seqs_to_align);

	int n_seqs_passing_swps3 = seqs_to_align->size();
	cout << "retained " << n_seqs_passing_swps3 << " seqs" << endl;

	if (justseqquery == true) { // if we are done, tell the user where the scores are written
		cout << "scores written to " << genefoldername + "/" + gene_name << ".seqquery" << endl;
		exit(0);
	}

	// TODO: assuming for now that all the user seqs are hits

	// get rid of duplicates
	remove_duplicates_SWPS3(seqs_to_align);
	int n_seqs_unique = seqs_to_align->size();
	cout << "removed " << n_seqs_passing_swps3 - n_seqs_unique << " taxonomic duplicates, leaving " << n_seqs_unique << " seqs to be aligned" << endl << endl;

	// if userguidetree overlaps less than a certain percentage, usertree = false
	if (usertree == true && userskipsearchdb == false) {
		double overlap = get_usertree_keepseq_overlap(seqs_to_align);
		overlap = 0.5;    // TODO: take this out - ? ceh -> make this user-settable?
		if (overlap < 0.5) {
			usertree = false;
			userguidetree = NULL;
			cerr << "user guide tree has insufficient overlap (" << overlap << ") and cannot be used. phlawd will quit now" << endl;
			exit(0);
		}
	}

	// get the list of files and record the higher taxa name and
	// add the additional sequences to the right hierarchy

	// we use two vectors here because for some purposes, we need only
	// one (such as when we have a user guide tree, we just need names),
	// while for others we need both (if we use the ncbi taxonomy we use
	// the ids and the names). see saturation_tests() for details on use

	vector<string> search_clade_tax_ids;
	vector<string> search_clade_tax_names;

	if (updateDB == true) {
		vector<int> toremove;
		for (int j = 0; j < seqs_to_align->size(); j++) {
			if (count(preexisting_db_seqs_ncbi_ids.begin(), preexisting_db_seqs_ncbi_ids.end(), seqs_to_align->at(j).get_ncbi_tax_id()) > 0) {
				toremove.push_back(j);
			}
		}
		for (int j = 0; j < toremove.size(); j++) {
			seqs_to_align->erase(seqs_to_align->begin() + toremove[toremove.size() - (j + 1)]);
		}
		//end the program if there is nothing to update
		if (seqs_to_align->size() == 0) {
			cout << "There are no new DB sequences to add." << endl;
			gifile.close();
		}
		cout << "total size of seqs to add (update):" << seqs_to_align->size() << endl;
		for (int j = 0; j < seqs_to_align->size(); j++) {
			cout << "adding: " << seqs_to_align->at(j).get_ncbi_tax_id() << endl;
		}
		cout << "adding sequences to db " << gene_db.get_db_name() << endl;
		gene_db.add_sequences(seqs_to_align);
		//add the right gi numbers before add the rest of the seqs are added to keep_seqs
		write_gi_numbers(seqs_to_align);
		gifile.close();
		//adding the update for user seqs
		toremove.clear();
		for (int j = 0; j < user_seqs->size(); j++) {
			if (count(preexisting_user_seqs_names.begin(), preexisting_user_seqs_names.end(), user_seqs->at(j).get_taxon_name()) > 0) {
				toremove.push_back(j);
			}
		}
		for (int j = 0; j < toremove.size(); j++) {
			user_seqs->erase(user_seqs->begin() + toremove[toremove.size() - (j + 1)]);
		}
		if (user_seqs->size() == 0) {
			cout << "There are no new user sequences to add." << endl;
			ufafile.close();
		}
		cout << "total size of user updated:" << user_seqs->size() << endl;
		for (int j = 0; j < user_seqs->size(); j++) {
			cout << "adding: " << user_seqs->at(j).get_taxon_name() << endl;
		}

		cout << "adding sequences to database: " << gene_db.get_db_name() << endl;
		gene_db.add_sequences(user_seqs);

		if (user_seqs->size() + seqs_to_align->size() == 0)
			exit(0);
		//if update is more than 20% then redo run
		if (user_seqs->size() + seqs_to_align->size() >= (preexisting_seqs.size() * MAX_UPDATE_PROPORTION)) {
			updateDB = false;
			cout << "update more than 20% of original size, redoing the whole thing" << endl;
			return 2;
		}
		//add the right gi numbers before add the rest of the seqs are added to keep_seqs
		write_user_numbers();
		ufafile.close();
		//check files for existing taxonomic break down as it will generally be
		//these seperations or more fine
		vector<string> align_names;
		gene_db.load_orig_alignment_labels_into(align_names); // in this case these are the ncbi ids
		if (usertree == false) {
			for (unsigned int i = 0; i < align_names.size(); i++) {
				string sql = "SELECT ncbi_id,left_value,right_value FROM taxonomy where ncbi_id = " + align_names[i] + ";";
				Query query(conn);
				query.get_result(sql);
				string t_id;
				string l_id;
				string r_id;
				while (query.fetch_row()) {
					t_id = query.getstr();
					l_id = query.getstr();
					r_id = query.getstr();
				}
				if (t_id.size() == 0) {
					cout << "There is an error getting the id for " << align_names[i] << endl;
					exit(0);
				}
				for (unsigned int j = 0; j < seqs_to_align->size(); j++) {
					//start here, need to get the higher taxon
					sql = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = '" + seqs_to_align->at(j).get_ncbi_tax_id() + "';";
					Query query2(conn);
					query2.get_result(sql);
					string lefttid;
					string righttid;
					while (query2.fetch_row()) {
						lefttid = query2.getstr();
						righttid = query2.getstr();
					}
					if (lefttid > l_id && righttid < r_id) {
						if ((int) count(search_clade_tax_ids.begin(), search_clade_tax_ids.end(), t_id) == 0) {
							search_clade_tax_ids.push_back(t_id);
							search_clade_tax_ids.push_back(align_names[i]);
							//add the sequences from the file into keep_seqs , this should be easier when moving to sqlite
							add_seqs_from_db_to_seqs_vector(align_names[i], seqs_to_align, preexisting_seqs);
							gene_db.remove_original_alignment_by_name(align_names[i]);
							cout << "deleting " << align_names[i] << endl;
						}
						break;
					}
				}
			}
			cout << "seqs to process: " << seqs_to_align->size() << endl;
		} else {	//usertree == true
					//can do seq similarity, tree building distance, or ncbi
					//going to try ncbi first
					//basically going to ask what are the ncbi taxa in each, if parent id of any of these is the same
					//going to test seq similarity and take best, if fails, then will just add to singletons

			//need to check if the tree is the same as the previous
			//if the nodes are then I can rename the new tree nodes to those
			vector<string> rem_files;
			//basic idea is to take the whole set of sequence names and see if the mrca of the names in a file include any seqs outside of that clade
			//storedseqs are the seqs from before
			FastaUtil fu;
			map<Node *, string> rename_nodes;	    //this should allow for renaming
			//the file_names left over are all the ones that need to be redone
			for (int j = 0; j < align_names.size(); j++) {
				bool test = true;
				vector<Sequence> tempseqs;
				bool is_aligned = true;
				string this_aln_fname = gene_name + "/" + align_names[j];
				fu.read_aligned_fasta_into(tempseqs, this_aln_fname);
				//TODO: if one seq need to just delete
				vector<string> mrca_names;
				//need to correct the names for user or not
				for (int i = 0; i < tempseqs.size(); i++) {
					mrca_names.push_back(tempseqs[i].get_taxon_name());
				}
				Node * tmrca = userguidetree->getMRCA(mrca_names);
				if (tmrca == NULL) {
					test = false;
				} else {
					vector<Node *> leaves = tmrca->get_leaves();
					for (int i = 0; i < leaves.size(); i++) {
						if ((int) count(mrca_names.begin(), mrca_names.end(), leaves[i]->getComment()) == 0
								&& (int) count(preexisting_db_seqs_ncbi_ids.begin(), preexisting_db_seqs_ncbi_ids.end(), leaves[i]->getComment()) > 0) {
							test = false;
							break;
						}
					}
				}
				if (test == true) {
					tmrca->setName(align_names[j]);
					rename_nodes[tmrca] = align_names[j];
				} else {
					add_seqs_from_db_to_seqs_vector(align_names[j], seqs_to_align, preexisting_seqs);
					gene_db.remove_original_alignment_by_name(align_names[j]);
					rem_files.push_back(align_names[j]);
				}
			}
			//remove
			for (int i = 0; i < rem_files.size(); i++) {
				vector<string>::iterator it;
				it = find(align_names.begin(), align_names.end(), rem_files[i]);
				//align_names.erase(it);
			}
			if (align_names.size() > 0) {            //explode is just a redo
				for (unsigned int i = 0; i < seqs_to_align->size(); i++) {
					int bestind;
					int bestscore = 0;
					for (unsigned int j = 0; j < align_names.size(); j++) {
						vector<Sequence> tempseqs;
						bool is_aligned = true;
						fu.read_aligned_fasta_into(tempseqs, align_names[j]);
						int tscore = get_single_to_group_seq_score(seqs_to_align->at(i), tempseqs);
						if (tscore > bestscore)
							bestind = j;
					}
					if ((int) count(search_clade_tax_names.begin(), search_clade_tax_names.end(), align_names[bestind]) == 0) {
						search_clade_tax_names.push_back(align_names[bestind]); //need to redo this file
						add_seqs_from_db_to_seqs_vector(align_names[bestind], seqs_to_align, preexisting_seqs);
						gene_db.remove_original_alignment_by_name(align_names[bestind]);
					}
				}
			} else {
				//rename the nodes with the remapped , the others can be random
				int ccount = 0;
				for (int i = 0; i < userguidetree->getInternalNodeCount(); i++) {
					if (rename_nodes.count(userguidetree->getInternalNode(i)) == 1) {
						userguidetree->getInternalNode(i)->setName(rename_nodes[userguidetree->getInternalNode(i)]);
					} else {
						while ((int) count(align_names.begin(), align_names.end(), to_string(ccount)) > 0) {
							ccount += 1;
						}
						userguidetree->getInternalNode(i)->setName(to_string(ccount));
					}
				}
				// setup for restart
				search_clade_tax_names.push_back(userguidetree->getRoot()->getName());
			}
		}
	} // end update prep

	// saturation tests
	cout << "starting saturation tests: identifying original clades to align before profiling" << endl;
	if (updateDB == true) {
		saturation_tests(search_clade_tax_ids, search_clade_tax_names, seqs_to_align);
	} else {
		write_gi_numbers(seqs_to_align);
		gifile.close();
		cout << "adding " << seqs_to_align->size() << " sequences to db " << gene_db.get_db_name() << endl;
		gene_db.add_sequences(seqs_to_align);
		if (userfasta == true) {
			write_user_numbers();
			ufafile.close();
			cout << "adding user sequences to database: " << gene_db.get_db_name() << endl;
			gene_db.add_sequences(user_seqs);
		}

		////////////////////////////////////////////////////////////////////////////////////////
		//
		//	tested to here in non-update case. sequence processing seems to be working.
		//  need to test update case.
		//
		////////////////////////////////////////////////////////////////////////////////////////

		search_clade_tax_ids.push_back(s_search_clade_id);
		search_clade_tax_names.push_back(clade_name);
		saturation_tests(search_clade_tax_ids, search_clade_tax_names, seqs_to_align);
	}

	logfile.close();
	delete known_seqs;
	delete seqs_to_align;
	return 0;
}

string SQLiteConstructor::get_cladename() {
	return clade_name;
}

vector<string> SQLiteConstructor::get_search() {
	return search;
}

string SQLiteConstructor::get_genename() {
	return gene_name;
}

double SQLiteConstructor::get_madcutoff() {
	return mad_cutoff;
}

double SQLiteConstructor::get_coverage() {
	return coverage;
}

double SQLiteConstructor::get_identity() {
	return identity;
}

int SQLiteConstructor::get_numthreads() {
	return numthreads;
}

//private

void SQLiteConstructor::load_info_for_all_seqs_matching_search_into(vector<seqRecordTuple> & found_seqs) {

	/* finds all sequences *from the source db* with a description matching anything
	 * in the member variable search, and returns a vector of objects of type
	 * seq_rec_tuple, which contain the sqlite row id, and taxon id of each
	 * sequence found.
	 */

	Database conn(db);
	string sql;
	if (searchliteral) {
		sql = "SELECT id, ncbi_id FROM sequence WHERE " + search[0] + " ;";
	} else {
		if (search.size() == 1) {
			sql = "SELECT id, ncbi_id FROM sequence WHERE description LIKE '%" + search[0] + "%'";
		} else {
			sql = "SELECT id, ncbi_id FROM sequence WHERE ";
			for (int i = 0; i < search.size() - 1; i++) {
				sql = sql + "description LIKE '%" + search[i] + "%' OR ";
			}
			sql = sql + "description LIKE '%" + search[search.size() - 1] + "%';";
		}
	}
	//string sql = "SELECT * FROM bioentry WHERE description LIKE '%"+search+"%'";
//	sql = "SELECT * FROM bioentry WHERE description LIKE '%"+search+"%' OR description LIKE '%trnK%'"
	Query query(conn);
	query.get_result(sql);
	int co = 1;
	while (query.fetch_row()) {

		seqRecordTuple seq;
		seq.db_record_id = to_string(query.getval());
		seq.ncbi_tax_id = to_string(query.getval());
		//		vector<string> vals;
//		vals.push_back(to_string());
//		vals.push_back();
		found_seqs.push_back(seq);
		co += 1;
	}
	cout << "search number " << co << endl;
	query.free_result();
}

vector<Sequence> SQLiteConstructor::make_seqs_from_seq_tuples_for_taxon(int taxon_ncbi_id, vector<seqRecordTuple> & seq_tuples) {

	/* converts each seqRecordTuple (these are structs that contain only the
	 * ncbi_tax_id and sqlite_id for a sequence) into actual sequence objects
	 * as long as they represent child taxa of the pass taxon_ncbi_id.
	 *
	 * first we look up the ncbi taxon for each sequence tuples, and then
	 * check to see if it represents a child of the taxon identified by
	 * taxon_ncbi_id. if it is, then the full information on this sequence
	 * is extracted from the db and stored in a Sequence object. These Sequence
	 * objects are retained within a vector that is returned after all the
	 * sequences are processed.
	 */

	vector<Sequence> matched_seqs;
	Database conn(db);
	Query query(conn);

	// get the left and right values for the parent taxon
	int parent_left_value, parent_right_value;
	string sql = "SELECT left_value, right_value FROM taxonomy WHERE ncbi_id = " + to_string(taxon_ncbi_id);
	query.get_result(sql);
	while (query.fetch_row()) {
		parent_left_value = query.getval();
		parent_right_value = query.getval();
		main_left = parent_left_value;
		main_right = parent_right_value;
	}
	query.free_result();

	// set up variables to hold info about each seq's taxon
	string this_seq_ncbi_tax_id;
	int this_tax_db_id;
	int this_tax_left_value;
	int this_tax_right_value;
	string this_tax_edited_name;

	for (int i = 0; i < seq_tuples.size(); i++) {

		// first look up this seq's taxon
		Query qtax(conn);
		this_seq_ncbi_tax_id = seq_tuples[i].ncbi_tax_id;
		sql = "SELECT id, left_value, right_value, edited_name FROM taxonomy WHERE ncbi_id == " + this_seq_ncbi_tax_id
				+ " and name_class == 'scientific name';";
		qtax.get_result(sql);
		while (qtax.fetch_row()) {
			this_tax_db_id = qtax.getval();
			this_tax_left_value = qtax.getval();
			this_tax_right_value = qtax.getval();
			this_tax_edited_name = qtax.getstr();
		}
		qtax.free_result();

		// if the seq represents a child of our parent taxon, get its info
		if (this_tax_left_value > parent_left_value && this_tax_right_value < parent_right_value) {

			string this_seq_db_id = seq_tuples[i].db_record_id;

			// set up variables for sequence info
			string this_seq_accession_number;
			string this_seq_gi_number;
			string this_seq_description;
			string this_seq_sequence_unaligned;

			// get the info from the db
			Query qseq(conn);
			sql = "SELECT accession_id, identifier, description, seq FROM sequence WHERE id == " + this_seq_db_id;
			qseq.get_result(sql);
			while (qseq.fetch_row()) {
				this_seq_accession_number = qseq.getstr();
				this_seq_gi_number = qseq.getstr();
				this_seq_description = qseq.getstr();
				this_seq_sequence_unaligned = qseq.getstr();
			}
			qseq.free_result();

			// make the new sequence object, set the dbid (as an int, not string)
			Sequence tseq = Sequence();
			tseq.set_sqlite_id(atoi(this_seq_db_id.c_str()));

			// set other fields from the db
//			tseq.set_label(this_seq_ncbi_tax_id);
			tseq.set_unaligned_sequence(this_seq_sequence_unaligned);
			tseq.set_ncbi_tax_id(this_seq_ncbi_tax_id);
			tseq.set_ncbi_gi_number(this_seq_gi_number);
			tseq.set_taxon_name(this_tax_edited_name);

			// add it to the vector
			matched_seqs.push_back(tseq);
		}
	}
	return matched_seqs;
}

/*
 * this takes the literal name from the file so this should be prefiltered
 * to be something that ncbi will match with the names not the edited names
 */
vector<Sequence> SQLiteConstructor::use_only_names_from_file(vector<Sequence> & seqs) {

	// TODO: clean this up...

	SQLiteDBController dbc = SQLiteDBController(db);
	Database conn(db);
	vector<string> * taxa = new vector<string>();
	vector<string> * taxa_ids = new vector<string>();
	//read file
	ifstream ifs(listfilename.c_str());
	string line;
	while (getline(ifs, line)) {
		TrimSpaces(line);
		taxa->push_back(line);
		string sql = "SELECT ncbi_id FROM taxonomy WHERE name = '" + line + "'";
		Query query1(conn);
		query1.get_result(sql);
		while (query1.fetch_row()) {
			int id;
			id = query1.getval();
			taxa_ids->push_back(to_string(id));
		}
		query1.free_result();
	}

	if (containswild)
		cout << "taxa (including all children) to be included: ";
	else
		cout << "exact taxon names to be matched (not including children): ";

	// just print a list of the names we're searching
	bool first = true;
	for (int i = 0; i < taxa_ids->size(); i++) {
		if (first)
			first = false;
		else
			cout << ", ";
		if (i > 30) {
			cout << "...";
			break;
		}
		int this_tax_id = atoi(taxa_ids->at(i).c_str());
		cout << dbc.get_sci_name_for_ncbi_tax_id(this_tax_id);
	}
	cout << endl;

	/*
	 * the adding of wild taxa (this will add ALL the children of the names in the list file)
	 */
	if (containswild == true) {
		vector<string> new_ids;
		for (int i = 0; i < taxa_ids->size(); i++) {
			//first get the left and right value for the taxa
			string sql = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = " + taxa_ids->at(i) + ";";
			Query query2(conn);
			query2.get_result(sql);
			string lefts;
			string rights;
			while (query2.fetch_row()) {
				int left = query2.getval();
				int right = query2.getval();
				lefts = to_string(left);
				rights = to_string(right);
			}
			sql = "SELECT ncbi_id FROM taxonomy WHERE left_value > " + lefts + " AND right_value < " + rights + " AND name_class = 'scientific name';";
			Query query(conn);
			query.get_result(sql);
			long count = query.num_rows();
			//exit(0);
			if (count == 0) {
				continue;
			} else {
				while (query.fetch_row()) {
					new_ids.push_back(to_string(query.getval()));
				}
			}
			query.free_result();
		}
		for (int i = 0; i < new_ids.size(); i++) {
			taxa_ids->push_back(new_ids[i]);
		}
	}
	/*
	 * end of the wild taxa
	 */

	cout << "sequences will be filtered against " << taxa_ids->size() << " taxon names" << endl;
	ifs.close();
	//end read file
	vector<Sequence> seqs_fn;
	for (int i = 0; i < seqs.size(); i++) {
		string taxid = seqs[i].get_ncbi_tax_id();
		int scount = count(taxa_ids->begin(), taxa_ids->end(), taxid);
		if (scount > 0) {
			seqs_fn.push_back(seqs[i]);
		}
	}
	/*
	 * added for higher taxa
	 */
	if (containshigher == true && containswild == false) {

		// TODO: what does this do? need to clarify output...

		cout << "this file contains higher taxa" << endl;
		for (int i = 0; i < taxa_ids->size(); i++) {
			string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = " + taxa_ids->at(i) + " and name_class = 'scientific name';";
			Query query(conn);
			query.get_result(sql);
			long count = query.num_rows();
			cout << count << endl;
			if (count == 0) {
				continue;
			} else {
				try {
					Sequence tse = find_best_exemplar_for_higher_taxon(taxa_ids->at(i), seqs);
					seqs_fn.push_back(tse);
				} catch (int a) {

				}
			}
			query.free_result();
		}
	}
	/*
	 * end added higher taxa
	 */
	delete taxa;
	delete taxa_ids;
	return seqs_fn;
}

/* just filters a set of sequences against a set of ncbi ids, saving the
 * ones corresponding to taxa within the set of ids into filtered_seqs.
 */
void SQLiteConstructor::load_sequences_filtered_by_ncbi_taxon_id_into(vector<Sequence> & filtered_seqs, vector<Sequence> & seqs_to_filter,
		vector<string> & ncbi_ids) {
	for (int i = 0; i < seqs_to_filter.size(); i++) {
		string taxid = seqs_to_filter[i].get_ncbi_tax_id();
		int scount = count(ncbi_ids.begin(), ncbi_ids.end(), taxid);
		if (scount > 0) {
			filtered_seqs.push_back(seqs_to_filter[i]);
		}
	}
}

/* called from use_only_names_from_file,  includes the procedure
 * to deal with higher taxa -- this does a mini phlawd construct
 * on the higher taxa to pick the best one
 * steps are
 * 1) send higher taxa name from use_only_names_from_file (which sends
 * based on a name in the list having children)
 * 2) get all seqs of the higher taxa from seqs (sent along)
 * 3) ortho check again known files
 * 4) take best seq (best overall)
 * 5) store in highertaxa container
 * 6) add these to those that pass the ortho check later (need to store
 * the seq in a file so that the higher taxa id is associated with the
 * smaller one (higher is in the phlawd file, smaller is store) saved
 * with the seq)
 */
//Sequence SQLiteConstructor::add_higher_taxa(string taxon_id, vector<Sequence> & all_seqs) {
Sequence SQLiteConstructor::find_best_exemplar_for_higher_taxon(string higher_taxon_id, vector<Sequence> & all_seqs) {

	/* ah, it seems that we are doing here is just finding the sequence that best
	 * represents the passed taxon, and adding it to the phlawd db under that taxon's
	 * name. this has the effect of creating a sequence in the phlawd db that is an
	 * exemplar for this higher taxon. */

	SQLiteDBController dbc = SQLiteDBController(db);
	string higher_taxon_name = dbc.get_sci_name_for_ncbi_tax_id(atoi(higher_taxon_id.c_str()));

	// find the seqs in the all_seqs vector that are within this higher taxon
	vector<string> child_tax_ids = get_valid_ncbi_child_taxon_ids_for_parent_id(higher_taxon_id);
	vector<Sequence> all_exemplar_seqs;
	load_sequences_filtered_by_ncbi_taxon_id_into(all_exemplar_seqs, all_seqs, child_tax_ids);

	if (all_exemplar_seqs.size() == 0) {
		cerr << "no sequences exist for ncbi taxon id " << higher_taxon_name << ". phlawd will exit";
		exit(0);
	}

	// apply the swps coverage/identity tests to make sure at least some representative seqs pass
	vector<Sequence> * potential_exemplar_seqs = new vector<Sequence>();
	get_best_hits_openmp_SWPS3(all_exemplar_seqs, potential_exemplar_seqs);

	if (potential_exemplar_seqs->size() < 0) {
		cerr << "all existing sequences for ncbi taxon id " << higher_taxon_name << " failed the coverage/identity tests. phlawd will exit";
		exit(0);
	}

	// score each seq in known_seqs against itself to get max scores
	vector<int> scores;
	SBMatrix mat = swps3_readSBMatrix("EDNAFULL");
	for (int i = 0; i < known_seqs->size(); i++) {
		scores.push_back(get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(i), &known_seqs->at(i)));
	}

	// find the exemplar that is the best match to something in known_seqs
	int bestid = 0;
	double bestiden = 0;
	for (int i = 0; i < potential_exemplar_seqs->size(); i++) {
		Sequence tseq = potential_exemplar_seqs->at(i);
		double maxiden = 0;
		bool rc = false;
		for (int j = 0; j < known_seqs->size(); j++) {
			bool trc = false;
			int ret = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &tseq);
			double tsc = double(ret) / double(scores[j]);

			/*					Sequence tseqrc = Sequence();
			 tseqrc.set_label(tseq.get_label()); */
			Sequence tseqrc = tseq.clone();

			tseqrc.set_unaligned_sequence(tseq.reverse_complement());
			int retrc = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &tseqrc);
			if (retrc > ret) {
				trc = true;
				tsc = double(retrc) / double(scores[j]);
			}
			if (tsc > maxiden) {
				maxiden = tsc;
			}
			//cout << tsc << endl;
		}
		if (maxiden >= bestiden) {
			bestid = i;
			bestiden = maxiden;
		}
	}
	Sequence bestseq = potential_exemplar_seqs->at(bestid);

	bestseq.set_ncbi_tax_id(higher_taxon_id);
	bestseq.set_taxon_name(higher_taxon_name);

	// TODO: clean out output
	cout << "higher taxa change" << endl;
	cout << higher_taxon_id << "=" << bestseq.get_ncbi_tax_id() << endl;
	logfile << "higher taxa change\n";
	logfile << "ncbi: " << bestseq.get_ncbi_tax_id() << "\n";
	logfile << "name: " << bestseq.get_taxon_name() << "\n";
	delete potential_exemplar_seqs;

	/*
	 * return the best sequence
	 */
	return bestseq;
}

/*
 * excluding taxa from sequences
 */
vector<Sequence> SQLiteConstructor::exclude_names_from_file(vector<Sequence>& seqs) {

	// TODO: clean this up...

	Database conn(db);
	vector<string> * taxa_ids = new vector<string>();
	//read file
	ifstream ifs(excludefilename.c_str());
	string line;
	vector<string> exclude_taxa;

	while (getline(ifs, line)) {
//		TrimSpaces(line); // don't want to trim spaces because we want to be able to isolate words as search terms
		exclude_taxa.push_back(line);
	}

	// just print a list of the names we're searching
	bool first = true;
	cout << "search terms for taxon names to be excluded (not including children): ";
	for (int i = 0; i < exclude_taxa.size(); i++) {
		if (first)
			first = false;
		else
			cout << ", ";

		if (i > 30) {
			cout << "...";
			break;
		}

		cout << "'" << exclude_taxa[i] << "'";
	}
	cout << endl;

	while (exclude_taxa.empty() == false) {
		string excl_search = exclude_taxa.back();
		exclude_taxa.pop_back();
		string sql;
		
		if (exclsearchliteral) { // accepts SQL syntax for searching for excluded names
		
            sql = "SELECT ncbi_id FROM taxonomy WHERE left_value > " + int_to_string(main_left) +
                    " AND right_value < " + int_to_string(main_right) +
                    " AND (" + excl_search + ") AND name_class == 'scientific name'";

        } else if (excl_search[0] == '*') { // uses wildcard searching to ignore any taxa with this in the name

            string excl_search_trimmed = excl_search.substr(1, excl_search.size());
            sql = "SELECT ncbi_id FROM taxonomy WHERE left_value > " + int_to_string(main_left) +
                    " AND right_value < " + int_to_string(main_right) +
                    " AND name like '%" + excl_search_trimmed + "%'  and name_class == 'scientific name'";

        } else {
            sql = "SELECT ncbi_id FROM taxonomy WHERE name = '" + excl_search + "'";

        }

		Query query1(conn);
		query1.get_result(sql);
		while (query1.fetch_row()) {
			int id;
			id = query1.getval();
			taxa_ids->push_back(to_string(id));
		}
		query1.free_result();
	}

	cout << "found " << taxa_ids->size() << " taxon names that will be excluded" << endl;
	ifs.close();
	//end read file
	vector<Sequence> seqs_fn;
	for (int i = 0; i < seqs.size(); i++) {
		string taxid = seqs[i].get_ncbi_tax_id();
		bool found = false;
		for (int j = 0; j < taxa_ids->size(); j++) {
			if (taxa_ids->at(j) == taxid) {
				found = true;
				break;
			}
		}
		if (found == false) {
			seqs_fn.push_back(seqs[i]);
		}
	}
	delete taxa_ids;
	return seqs_fn;
}

/*
 * excluding gis from sequences
 */
vector<Sequence> SQLiteConstructor::exclude_gis_from_file(vector<Sequence> &seqs) {
	vector<string> * gi_ids = new vector<string>();
	//read file
	ifstream ifs(exclude_gi_filename.c_str());
	string line;
	while (getline(ifs, line)) {
		TrimSpaces(line);
		gi_ids->push_back(line);
	}
	cout << gi_ids->size() << " gis in the file" << endl;
	ifs.close();
	//end read file
	vector<Sequence> seqs_fn;
	for (int i = 0; i < seqs.size(); i++) {
		string giid = seqs[i].get_ncbi_gi_number();
		bool found = false;
		for (int j = 0; j < gi_ids->size(); j++) {
			if (gi_ids->at(j) == giid) {
				found = true;
				break;
			}
		}
		if (found == false) {
			seqs_fn.push_back(seqs[i]);
		}
	}
	delete gi_ids;
	return seqs_fn;
}

/*
 * include only gis from file
 */
vector<Sequence> SQLiteConstructor::include_gis_from_file(vector<Sequence> & seqs) {
	vector<string> * gi_ids = new vector<string>();
	//read file
	ifstream ifs(include_gi_filename.c_str());
	string line;
	while (getline(ifs, line)) {
		TrimSpaces(line);
		gi_ids->push_back(line);
	}
	cout << gi_ids->size() << " gis in the file" << endl;
	ifs.close();
	//end read file
	vector<Sequence> seqs_fn;
	for (int i = 0; i < seqs.size(); i++) {
		//string giid = seqs[i].get_accession();
		string giid = seqs[i].get_ncbi_gi_number();
		bool found = false;
		for (int j = 0; j < gi_ids->size(); j++) {
			if (gi_ids->at(j) == giid) {
				found = true;
				break;
			}
		}
		if (found == true) {
			seqs_fn.push_back(seqs[i]);
		}
	}
	delete gi_ids;
	return seqs_fn;
}

/*
 * OPENMP version
 */
void SQLiteConstructor::get_best_hits_openmp_SWPS3_justquery(vector<Sequence> & seqs_to_score, vector<Sequence> * seqs_to_keep) {

	/* scores all the the seqs in seqs_to_score against known_seqs, and adds the ones
	 * scoring high enough (against any seq in known_seqs) to seqs_to_keep. */

	cout << "Starting cover/ident calculations" << endl;
	vector<int> known_scores;
	SBMatrix mat = swps3_readSBMatrix("EDNAFULL");
	//SBMatrix mat = swps3_get_premade_SBMatrix( "EDNAFULL" );

	for (int i = 0; i < known_seqs->size(); i++) {
		known_scores.push_back(get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(i), &known_seqs->at(i)));
	}
	map<Sequence*, bool> keep_seqs_rc_map;
	vector<double> justqueryvec;
	vector<double> justqueryvec2;
	vector<string> justqueryname;

	omp_set_dynamic(0);
	omp_set_num_threads(get_numthreads());
	#pragma omp parallel for shared(keep_seqs_rc_map,justqueryvec,justqueryvec2,justqueryname)
	for (int i = 0; i < seqs_to_score.size(); i++) {
		double maxide = 0;
		double maxcov = 0;
		bool rc = false;
		for (int j = 0; j < known_seqs->size(); j++) {
//        cout << seqs[i].get_id() << endl;
			bool trc = false;
			int ret = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &seqs_to_score[i]);
			double tsc = double(ret) / double(known_scores[j]);
			seqs_to_score[i].perm_reverse_complement(); //make reverse complement
			int retrc = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &seqs_to_score[i]);
			seqs_to_score[i].perm_reverse_complement(); //make it back to original
			int setsc = get_swps3_score_and_rc_cstyle(mat, &seqs_to_score[i], &seqs_to_score[i]);
			double fsetsc = double(ret) / double(setsc);
			if (retrc > ret) {
				trc = true;
				tsc = double(retrc) / double(known_scores[j]);
				fsetsc = double(retrc) / double(setsc);
			}
//	    cout <<i << " " << j << " " << setsc << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << " " << double(ret)/double(setsc) << endl;
//	    cout << seqs[i].get_sequence() << endl;
//	    cout << known_seqs->at(j).get_sequence() << endl;
//	    cout << seqs[i].get_sequence().size() << endl;
//	    cout << known_seqs->at(j).get_sequence().size() << endl;
//	    exit(0);
			if (tsc > maxide && std::numeric_limits<double>::infinity() != tsc) {
				maxide = tsc;
				maxcov = fsetsc;
				rc = trc;
			}
		}
		justqueryvec.push_back(maxide);
		justqueryvec2.push_back(maxcov);
		justqueryname.push_back(seqs_to_score[i].get_ncbi_gi_number());
		if (maxide >= identity && maxcov >= coverage) {
			keep_seqs_rc_map[&seqs_to_score[i]] = rc;
			//with this we don't have to keep track of rc anymore unless we want to
			if (rc == true)
				seqs_to_score[i].perm_reverse_complement(); //the sequence is suppose to be reverse complement
		}
	}

	map<Sequence*, bool>::iterator it;
	for (it = keep_seqs_rc_map.begin(); it != keep_seqs_rc_map.end(); it++) {
		seqs_to_keep->push_back(*(*it).first);
	}
	ofstream outfile;
	outfile.open((gene_name + ".seqquery").c_str(), ios::out);
	for (int i = 0; i < justqueryvec.size(); i++) {
		outfile << justqueryname[i] << "\t" << justqueryvec[i] << "\t" << justqueryvec2[i] << endl;
	}
	outfile.close();
}

void SQLiteConstructor::get_best_hits_openmp_SWPS3(vector<Sequence> & seqs, vector<Sequence> * seqs_to_keep) {
	vector<int> known_scores;
	SBMatrix mat = swps3_readSBMatrix("EDNAFULL");
	//SBMatrix mat = swps3_get_premade_SBMatrix( "EDNAFULL" );
	for (int i = 0; i < known_seqs->size(); i++) {
		known_scores.push_back(get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(i), &known_seqs->at(i)));
	}
	map<Sequence*, bool> keep_seqs_rc_map;
	
	omp_set_dynamic(0);
	omp_set_num_threads(get_numthreads());
	#pragma omp parallel for shared(keep_seqs_rc_map)
	for (int i = 0; i < seqs.size(); i++) {
		double maxide = 0;
		double maxcov = 0;
		bool rc = false;
		for (int j = 0; j < known_seqs->size(); j++) {
			bool trc = false;
			int ret = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &seqs[i]);
			double tsc = double(ret) / double(known_scores[j]);
			seqs[i].perm_reverse_complement(); //make reverse complement
			int retrc = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &seqs[i]);
			seqs[i].perm_reverse_complement(); //make it back to original
			int setsc = get_swps3_score_and_rc_cstyle(mat, &seqs[i], &seqs[i]);
			double fsetsc = double(ret) / double(setsc);
//	    cout <<i << " " << j << " " << setsc << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << " " << double(ret)/double(setsc) << endl;
//	    cout << seqs[i].get_sequence() << endl;
//	    cout << known_seqs->at(j).get_sequence() << endl;
//	    cout << seqs[i].get_sequence().size() << endl;
//	    cout << known_seqs->at(j).get_sequence().size() << endl;
//	    exit(0);
			if (retrc > ret) {
				trc = true;
				tsc = double(retrc) / double(known_scores[j]);
				fsetsc = double(retrc) / double(setsc);
			}
			if (tsc > maxide && std::numeric_limits<double>::infinity() != tsc) {
				maxide = tsc;
				maxcov = fsetsc;
				rc = trc;
			}
			if (maxide >= min(identity + (identity * 0.5), .99))
				break;
		}
		if (maxide >= identity && maxcov >= coverage) {
			keep_seqs_rc_map[&seqs[i]] = rc;
			//with this we don't have to keep track of rc anymore unless we want to
			if (rc == true)
				seqs[i].perm_reverse_complement(); //the sequence is suppose to be reverse complement
		}
	}
	
	map<Sequence*, bool>::iterator it;
	for (it = keep_seqs_rc_map.begin(); it != keep_seqs_rc_map.end(); it++) {
		seqs_to_keep->push_back(*(*it).first);
	}
}

/* adds guide sequences with taxon names found in the gb taxonomy to the set of sequence
 * to be aligned.
 * 
 * this is wildly inefficient in its current implementation. should be refactored to
 * use hashset data structures for fast lookups. */
/*void SQLiteConstructor::replace_gb_seqs_with_guide_seqs(vector<Sequence> * known_seqs, vector<Sequence> * in_seqs) {
    
    // get names of all guide seqs
    for (unsigned int i = 0; i < known_seqs->size(); i++) {
        
        in_seqs.push_back(known_seqs->at(i));
    }
    
/*    // compare gb seqs to guide seqs
	for (unsigned int i = in_seqs->size() - 1; i >= 0 ; i--) {

        string cur_name = in_seqs->at(i).get_taxon_name();

        // first check if we've seen this name yet
        if (known_seqs_names_saved.find(cur_name) != known_seqs_names.end()) {

            // we have seen this name in the guide seqs, erase this gb seq
            in_seqs->erase(in_seqs->begin() + i - 1);

        } else {
        
            // have not seen name yet, check if it is in the guide seqs
            if (known_seqs_names.find(cur_name) != known_seqs_names.end()) {
                
                // name is in guide seqs, erase this gb seq
                known_seqs_names_saved.insert(cur_name);
                in_seqs->erase(in_seqs->begin() + i - 1);

            }            
        }
    }

    for (unsigned int i = 0; i < known_seqs->size(); i++) {
        in_seqs->push_back(known_seqs->at(i));
    } 
} */

void SQLiteConstructor::remove_duplicates_SWPS3(vector<Sequence> * keep_seqs) {
	vector<string> ids;
	vector<string> unique_ids;
	int mycount;

	//uses NCBI taxon ids for dups

	for (unsigned int i = 0; i < keep_seqs->size(); i++) {
		ids.push_back(keep_seqs->at(i).get_ncbi_tax_id());
		mycount = 0;
		if (unique_ids.size() > 0) {
			mycount = (int) count(unique_ids.begin(), unique_ids.end(), keep_seqs->at(i).get_ncbi_tax_id());
		}
		if (mycount == 0) {
			unique_ids.push_back(keep_seqs->at(i).get_ncbi_tax_id());
		}
	}

	/*
	 * get the best score for each known seq
	 */
	vector<int> scores;
	SBMatrix mat = swps3_readSBMatrix("EDNAFULL");
	for (int i = 0; i < known_seqs->size(); i++) {
		scores.push_back(get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(i), &known_seqs->at(i)));
	}

	vector<int> remove;
	for (unsigned int i = 0; i < unique_ids.size(); i++) {
		int mycount = 0;
		mycount = (int) count(ids.begin(), ids.end(), unique_ids[i]);
		if (mycount > 1) {
			vector<int> tremove;
			for (int j = 0; j < ids.size(); j++) {
				if (ids[j] == unique_ids[i]) {
					remove.push_back(j);
					tremove.push_back(j);
				}
			}
			int bestid = 0;
			double bestiden = 0;
			for (int j = 0; j < tremove.size(); j++) {
				Sequence tseq = keep_seqs->at(tremove[j]);
				double maxiden = 0;
				for (int j = 0; j < known_seqs->size(); j++) {
					int ret = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &tseq);
					double tsc = double(ret) / double(scores[j]);
					if (tsc > maxiden) {
						maxiden = tsc;
					}
				}
				if (maxiden >= bestiden) {
					bestid = tremove[j];
					bestiden = maxiden;
				}

			}
			vector<int>::iterator it;
			it = find(remove.begin(), remove.end(), bestid);
			remove.erase(it);
		}
	}
	sort(remove.begin(), remove.end());
	reverse(remove.begin(), remove.end());
	for (unsigned int i = 0; i < remove.size(); i++) {
		keep_seqs->erase(keep_seqs->begin() + remove[i]);
	}
}

void SQLiteConstructor::reduce_genomes(vector<Sequence> * keep_seqs) {
	/*
	 * get the best score for each known seq
	 */
	vector<int> scores;
	SBMatrix mat = swps3_readSBMatrix("EDNAFULL");
	for (int j = 0; j < known_seqs->size(); j++) {
		scores.push_back(get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &known_seqs->at(j)));
	}
	for (unsigned int i = 0; i < keep_seqs->size(); i++) {
		if (keep_seqs->at(i).get_sequence().size() > shrinkablethreshold) { // TODO: make this user-settable
			cout << "sequence for " << keep_seqs->at(i).get_taxon_name() << " (gi " << keep_seqs->at(i).get_ncbi_gi_number() << ") " << " is longer than "
					<< shrinkablethreshold << ". it will be shrunk " << endl;
			Sequence tseq = keep_seqs->at(i);
			double maxiden = 0;
			int maxknown = 0;
			for (int j = 0; j < known_seqs->size(); j++) {
				int ret = get_swps3_score_and_rc_cstyle(mat, &known_seqs->at(j), &tseq);
				double tsc = double(ret) / double(scores[j]);
				if (tsc > maxiden) {
					maxiden = tsc;
					maxknown = j;
				}
			}
			const string tempfile = genefoldername + "genome_shrink_in";
			vector<Sequence> sc1;
			FastaUtil seqwriter;
			sc1.push_back(keep_seqs->at(i));
			//for (int j=0;j<known_seqs->size();j++){
			sc1.push_back(known_seqs->at(maxknown));
			//}
			seqwriter.writeFileFromVector(tempfile, sc1);
			string cmd = "mafft ";
			cmd += genefoldername + "genome_shrink_in > ";
			cmd += genefoldername + "genome_shrink_aln 2> genome_shrink.mafftlog";
//			cout << cmd << endl;
			cout << "using mafft to align the genome to a known sequence" << endl;
			system(cmd.c_str());
			/*string cmd2 = "phyutility -clean 0.5 -in ";
			 cmd2 += genefoldername+"genome_shrink_aln -out ";
			 cmd2 += genefoldername+"genome_shrink_out";
			 cout << cmd << endl;
			 cout << "cleaning" << endl;
			 system(cmd2.c_str());
			 */
			//instead of phyutility
			cout << "excising aligned sites (retaining these)" << endl;
			clean_shrunken_genomes();
			/*
			 * reading in the sequencing and replacing
			 */
			FastaUtil seqreader;
			vector<Sequence> trimmed_genomes;
			string trimmed_genomes_output = genefoldername + "genome_shrink_out";

			seqreader.read_aligned_fasta_into(trimmed_genomes, trimmed_genomes_output);
			int n_trimmed_genomes = trimmed_genomes.size();

			for (int j = 0; j < n_trimmed_genomes; j++) {
				if (trimmed_genomes.at(j).get_taxon_name() == keep_seqs->at(i).get_taxon_name()) {

					//replace "-" with ""
					string cleanSeq = trimmed_genomes.at(j).get_sequence();
					replaceAll(cleanSeq, "-", "");

					keep_seqs->at(i).set_aligned_sequence(cleanSeq);
				}
			}
			cout << "reduced size: " << keep_seqs->at(i).get_aligned_length() << endl << endl;
		}
	}
}

void SQLiteConstructor::clean_shrunken_genomes() {
	double threshold = 0.25; // missing more than this, then remove // this was set to 0.5 but then no sites are removed
	FastaUtil fu;
	vector<Sequence> tempalseqs;
	string trimmed_genomes = genefoldername + "genome_shrink_aln";
	fu.read_aligned_fasta_into(tempalseqs, trimmed_genomes); // assuming this is aligned

	size_t seqlength = tempalseqs[0].get_sequence().size();
	float fseql = float(tempalseqs.size());
	vector<size_t> removeem;
	for (size_t j = 0; j < seqlength; j++) {
		size_t gaps = 0;
		for (size_t i = 0; i < tempalseqs.size(); i++) {
			if (tempalseqs[i].get_sequence()[j] == '-' || tempalseqs[i].get_sequence()[j] == 'N' || tempalseqs[i].get_sequence()[j] == 'n')
				gaps += 1;
		}
		double curp = gaps / fseql;
		if (curp > threshold) {
			removeem.push_back(j);
		}
	}
	for (size_t i = 0; i < tempalseqs.size(); i++) {
		string a;
		for (size_t j = 0; j < seqlength; j++) {
			if (count(removeem.begin(), removeem.end(), j) == 0)
				a += tempalseqs[i].get_sequence()[j];
		}
		tempalseqs[i].set_aligned_sequence(a);
	}
	remove((genefoldername + "genome_shrink_aln").c_str());
	fu.writeFileFromVector(genefoldername + "genome_shrink_out", tempalseqs);
}

vector<string> SQLiteConstructor::get_valid_ncbi_child_taxon_ids_for_parent_id(string parent_taxon_id) {

	/* recursively walk the ncbi taxonomy down from the taxon identified
	 * the passed ncbi_taxon_id, using a vector to keep track of taxon ids
	 * we haven't yet checked. we start by adding the passed parent_taxon_id
	 * to the traversal vector, and add every child we find to both that
	 * vector, as well as a vector that remembers all the taxon ids we've
	 * seen. then we filter the latter vector to remove environmental samples,
	 * and return the vector containing the filtered set
	 *
	 * note: we do include environmental samples in the traversal, but then
	 * in the next step we exclude all the environmental samples we have found.
	 * it is not clear why they are included in the first place. one possible
	 * reason to include them in the traversal is that 'environmental samples'
	 * taxa can have children, and if their children might not be environmental
	 * samples, we would need to include those parents in order to find those
	 * children. not sure if this is true though.
	 */

	if (parent_taxon_id == "1") {
		cout << "special case: specified taxon (ncbi id = 1) is root for entire tree of life; getting ALL taxa" << endl;

		Database conn(db);
		vector<string> allids;
		string sql = "SELECT ncbi_id FROM taxonomy where name_class == 'scientific name' and ";
		sql += "name_class not like '%environmental%' and name not like '%environmental%';";
		Query query(conn);
		query.get_result(sql);

		while (query.fetch_row()) {

//			string this_name = query.getstr();
//			string this_name_class = query.getstr();
			string this_ncbi_taxon_id = query.getstr();

//			bool is_scientific_name = this_name_class.find("scientific") != string::npos;
//			bool is_not_environmental_seq = !(this_name.find("environmental") == string::npos || this_name_class.find("environmental") == string::npos);

//			if (is_scientific_name && is_not_environmental_seq) {
			allids.push_back(this_ncbi_taxon_id);
//			}
		}
		query.free_result();
		return allids;

	} else { // standard case; this is not the root node

		vector<string> all_child_taxon_ids;			// will hold the ncbi taxon ids of all taxa that we find, including the deepest parent itself
		vector<string> taxon_ids_to_check;			// will hold ncbi taxon ids that we haven't tried to find children for yet
		taxon_ids_to_check.push_back(parent_taxon_id);
		bool deepest_parent_is_tip = true;			// assuming the original parent_taxon_id is has no children until we find out otherwise

		Database conn(db);

		while (!taxon_ids_to_check.empty()) {

			// get a putative parent taxon id from the stack
			string this_parent_taxon_id = taxon_ids_to_check.back();
			taxon_ids_to_check.pop_back();

			// look for its children
			string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id == " + this_parent_taxon_id + " and name_class == 'scientific name';";
			Query query(conn);
			query.get_result(sql);

			// if we ever find any children, then the original parent must not be a tip
			if (query.num_rows() > 0) {
				deepest_parent_is_tip = false;
			}

			while (query.fetch_row()) {
				string this_ncbi_taxon_id = query.getstr();
				taxon_ids_to_check.push_back(this_ncbi_taxon_id);	// add this child to the list of taxa to check for further children
				all_child_taxon_ids.push_back(this_ncbi_taxon_id);	// remember that this is a child of the original parent
			}
			query.free_result();
		}

		// if there were no children; just save the passed parent_taxon_id itself
		if (deepest_parent_is_tip) {
			all_child_taxon_ids.push_back(parent_taxon_id);
		}

		// now look through the taxa we found and exclude all the environmental samples
		vector<string> validated_child_taxon_ids;
		// vector<string> validated_child_taxon_names; // this is never used
		while (!all_child_taxon_ids.empty()) {
			string this_child_taxon_id = all_child_taxon_ids.back();
			all_child_taxon_ids.pop_back();

			// get name information about this taxon
			string sql = "SELECT name FROM taxonomy WHERE ncbi_id == " + this_child_taxon_id + " and name_class == 'scientific name';";
			Query query(conn);
			query.get_result(sql);

			// save taxa that are not environmental samples
			while (query.fetch_row()) {
				string this_child_taxon_name = query.getstr();
				bool is_not_environmental_sample = (int) this_child_taxon_name.find("environmental") == string::npos ? true : false;

				if (is_not_environmental_sample) {
					validated_child_taxon_ids.push_back(this_child_taxon_id);

					// validated_child_taxon_names.push_back(this_child_taxon_name); // this is never used
				}
			}
			query.free_result();
		}

		/*		while (!taxon_ids_to_check.empty()) {
		 // attempt to get all children for this taxon
		 string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = " + taxon_ids_to_check.back() + " and name_class='scientific name';";
		 Query query(conn);
		 query.get_result(sql);
		 taxon_ids_to_check.pop_back();

		 if (query.num_rows() > 0) {
		 deepest_parent_is_tip = false;
		 }
		 while (query.fetch_row()) {

		 string this_ncbi_taxon_id = query.getstr();
		 taxon_ids_to_check.push_back(this_ncbi_taxon_id);
		 all_child_taxon_ids.push_back(this_ncbi_taxon_id);
		 }
		 query.free_result();
		 }



		 vector<string> validated_child_taxon_ids;
		 vector<string> allnames;

		 for (int i = 0; i < all_child_taxon_ids.size(); i++) {
		 string sql = "SELECT name,name_class FROM taxonomy WHERE ncbi_id = ";
		 sql += all_child_taxon_ids[i];
		 Query query(conn);
		 query.get_result(sql);
		 //StoreQueryResult R = query.store();
		 while (query.fetch_row()) {
		 //string tid = R[j][0].c_str();
		 string tn = query.getstr();
		 string cln = query.getstr();
		 if (cln.find("scientific") != string::npos && tn.find("environmental") == string::npos && cln.find("environmental") == string::npos) {
		 validated_child_taxon_ids.push_back(all_child_taxon_ids[i]); //was taxon id, now ncbi id
		 allnames.push_back(tn);
		 }
		 }
		 query.free_result();
		 } */
		return validated_child_taxon_ids; // TODO: would be better to use pass by reference
	}
}

/*
 * getting the final ncbi_ids for the children originating at a node
 */
vector<string> SQLiteConstructor::get_final_children_node(Node * node) {
	vector<string> allids;
	vector<Node *> leaves = node->get_leaves();
	for (int i = 0; i < leaves.size(); i++) {
		allids.push_back(leaves[i]->getComment());
	}
	return allids;
}

/*
 * same as get_final_children_node but added the ability to do hierarchical (so the tips in tree may or may not be species)
 */
vector<string> SQLiteConstructor::get_final_children_node_hier(Node * node) {
	vector<string> allids;
	vector<Node *> leaves = node->get_leaves();
	for (int i = 0; i < leaves.size(); i++) {
		vector<string> tempids = get_valid_ncbi_child_taxon_ids_for_parent_id(leaves[i]->getComment()); // TODO: why are we using the comment field?
		for (int j = 0; j < tempids.size(); j++) {
			allids.push_back(tempids[j]);
		};
	}
	return allids;
}

//void SQLiteConstructor::find_db_child_seqs_of_ncbi_taxon_id(string parent_id, vector<Sequence> * seqs_to_search, vector<Sequence> * found_seqs) {
	
void SQLiteConstructor::get_child_seqs_for_ncbi_tax_id(string parent_id, vector<Sequence> * seqs_to_search, vector<Sequence> * found_seqs) { //, vector<Sequence> * found_seqs_user,) {

	/* this function filters sequences in the seqs_to_filter vector against
	 * a list of valid child taxon ids for some parent provided by the function
	 * get_valid_ncbi_child_taxon_ids_for_ncbi_parent_taxon_id().
	 *
	 * seqs that pass the filtering step (i.e. the taxa they represent are within
	 * the list of valid child taxon ids for the taxon identified by
	 * parent_ncbi_id) are added to the vector filtered_seqs.
	 */

	// get the ncbi taxon ids of all the valid terminal children of this parent taxon
	vector<string> child_tax_ids = get_valid_ncbi_child_taxon_ids_for_parent_id(parent_id);
	set<string> user_seq_tax_ids; // remember the ids of user seqs we added so we don't add duplicates from the db

	// first get all the user sequences for this taxon
	for (unsigned int i = 0; i < user_seqs->size(); i++) {

		// look for this sequence's taxon id in the list of valid ids
		string user_seq_tax_id = user_seqs->at(i).get_ncbi_tax_id();
		bool this_seq_represents_a_child = (int) count(child_tax_ids.begin(), child_tax_ids.end(), user_seq_tax_id) > 0 ? true : false;

		if (this_seq_represents_a_child) {
			found_seqs->push_back(user_seqs->at(i));
			user_seq_tax_ids.insert(user_seq_tax_id);
		}
	}
	
	SQLiteDBController dbc = SQLiteDBController(db);

	// now add db seqs
	for (unsigned int i = 0; i < seqs_to_search->size(); i++) {

		// look for this sequence's taxon id in the list of valid ids
		string cur_seq_tax_id = seqs_to_search->at(i).get_ncbi_tax_id();
		bool this_seq_represents_a_child = (int) count(child_tax_ids.begin(), child_tax_ids.end(), cur_seq_tax_id) > 0 ? true : false;

		if (this_seq_represents_a_child) {
			
			// don't add the seq if its taxon is already represented by a user sequence
			if (user_seq_tax_ids.find(cur_seq_tax_id) != user_seq_tax_ids.end()) {
				cout << "excluding a db sequence for \"" << dbc.get_sci_name_for_ncbi_tax_id(atoi(cur_seq_tax_id.c_str())) <<
					"\" because it is already represented by a user sequence." << endl;

			} else { // no user sequence for this child taxon, so add it
				found_seqs->push_back(seqs_to_search->at(i));
			}
		}
	}
	
	// report to user
	cout << "found " << found_seqs->size() - user_seq_tax_ids.size() << " database sequences and " <<
		user_seq_tax_ids.size() << " user sequences for taxon \"" << dbc.get_sci_name_for_ncbi_tax_id(atoi(parent_id.c_str())) <<
			"\" (ncbi id = " << parent_id << ")" << endl;
}

/*
void SQLiteConstructor::get_seqs_for_names_user(string inname_id, vector<Sequence> * temp_seqs) {

	vector<string> final_ids;
	final_ids = get_valid_ncbi_child_taxon_ids_for_parent_id(inname_id);


} */

/* for userguidetree
 * this is intended to retrieve all the seqs that are contained below a node
 */
void SQLiteConstructor::get_seqs_for_nodes(Node * node, vector<Sequence> * seqs, vector<Sequence> * temp_seqs) {
	vector<string> final_ids;
	//final_ids = get_final_children_node(node);//TODO: see below
	final_ids = get_final_children_node_hier(node); //TODO: decide between these two
	for (unsigned int i = 0; i < seqs->size(); i++) {
		string tid = seqs->at(i).get_ncbi_tax_id();
		int mycount = 0;
		mycount = (int) count(final_ids.begin(), final_ids.end(), tid);
		if (mycount > 0) {
			temp_seqs->push_back(seqs->at(i));
		}
	}
}

void SQLiteConstructor::get_seqs_for_user_nodes(Node * node, vector<Sequence> * temp_seqs) {
	vector<string> final_ids;
	vector<Node *> leaves = node->get_leaves();
	for (unsigned int i = 0; i < user_seqs->size(); i++) {
		if (user_fasta_node_map.count(&user_seqs->at(i)) > 0) {
			int mycount = 0;
			mycount = (int) count(leaves.begin(), leaves.end(), user_fasta_node_map[&user_seqs->at(i)]);
			if (mycount > 0)
				temp_seqs->push_back(user_seqs->at(i));
		}
	}
}

void SQLiteConstructor::make_mafft_multiple_alignment(vector<Sequence> * inseqs) {

	/* just a wrapper for the version of this function that
	 * accepts two vectors; if we only have one then we just
	 * pass an empty vector for the second one.
	 *
	 * TODO: make the original function accept an arbitrary
	 * number of sequence vectors so we can just have one
	 * function that can combine things. or have another
	 * function that combines sequence vectors...
	 */

	vector<Sequence> emptys;
	make_mafft_multiple_alignment(inseqs, &emptys);
}

void SQLiteConstructor::make_mafft_multiple_alignment(vector<Sequence> * inseqs1, vector<Sequence> * inseqs2) {

	/* accepts two vectors of sequences (typically one for the db seqs
	 * and one for the user seqs), concatenates the vectors, and calls
	 * mafft to make an alignment for all the seqs.
	 */

	//make file
	vector<double> retvalues;
	FastaUtil seqwriter1;
	vector<Sequence> sc1;
	for (unsigned int i = 0; i < inseqs1->size(); i++) {
		sc1.push_back(inseqs1->at(i));
	}
	//TODO: do the user sequence names need to be changed
	for (unsigned int i = 0; i < inseqs2->size(); i++) {
		sc1.push_back(inseqs2->at(i));
	}
	const string fn1 = genefoldername + "tempfile";
	seqwriter1.writeFileFromVector(fn1, sc1);

	//make alignment
	string cmd = "mafft "; //--thread " + to_string(omp_get_max_threads());
	cmd += genefoldername + "tempfile > ";
	cmd += genefoldername + "outfile 2> ";
	cmd += genefoldername + "mafft.out";
//	cout << "aligning" << endl;
//	cout << cmd << endl;
	/*
	 FILE *fp = popen(cmd.c_str(), "r" );
	 char buff[1000];
	 while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	 string line(buff);
	 }
	 pclose( fp );
	 */
	system(cmd.c_str());
}

double SQLiteConstructor::calculate_MAD_quicktree() {

	/* use quicktree to calculate the MAD score. we feed it the alignment
	 * from mafft (assuming here that such an alignment exists) and it
	 * calculates pairwise taxon-taxon scores for all the taxa in the
	 * alignment. the MAD is calculated using this matrix. */

	string line;
	ofstream mafft_out_clean;

	// clean the mafft alignment and write back to 'mafft_out.clean'
	mafft_out_clean.open((genefoldername + "mafft_out.clean").c_str(), ios::out);
	FastaUtil fu;
	vector<Sequence> tempseqs;
	string mafft_output = genefoldername + "outfile";
	fu.read_aligned_fasta_into(tempseqs, mafft_output);
	for (int i = 0; i < tempseqs.size(); i++) {
		mafft_out_clean << tempseqs[i].get_taxon_name() << "\t" << tempseqs[i].get_aligned_seq() << endl;
	}
	mafft_out_clean.close();

	// use quicktree to make a sequence similarity matrix of all tips this alignment
	string cmd = "quicktree -in a -out m ";
	cmd += genefoldername + "mafft_out.clean > ";
	cmd += genefoldername + "quicktree_dists";
	cout << "calculating alignment pairwise distances" << endl;
	system(cmd.c_str());

	vector<double> p_values;
	vector<double> jc_values;

	/*
	 * read the matrix
	 */
	//string line;
	ifstream pfile((genefoldername + "quicktree_dists").c_str());
	vector<string> tokens;
	int nspecies = 0;
	int curspecies = 0;
	bool begin = true;
	if (pfile.is_open()) {
		while (!pfile.eof()) {
			getline(pfile, line);
			string del("\t ");
			tokens.clear();
			Tokenize(line, tokens, del);
			if (tokens.size() > 1) {
				double n1;
				for (int j = curspecies; j < nspecies - 1; j++) {
					n1 = atof(tokens.at(j + 2).c_str());
					p_values.push_back(n1);
					//jc will be (-(3./4.)*math.log(1-(4./3.)*p))
					jc_values.push_back((-(3. / 4.) * log(1 - (4. / 3.) * n1)));
				}
				curspecies += 1;
			} else if (begin == true) {
				begin = false;
				TrimSpaces(line);
				nspecies = atoi(line.c_str());
			}
		}
		pfile.close();
	}
	vector<double> all_abs_devs;
	double alignment_median_dev = 0;
	for (unsigned int i = 0; i < p_values.size(); i++) {
		all_abs_devs.push_back(fabs(jc_values[i] - p_values[i]));
	}
	alignment_median_dev = median(all_abs_devs);
	cout << "median deviation: " << alignment_median_dev << endl;
	vector<double> all_med_devs;
	for (unsigned int i = 0; i < p_values.size(); i++) {
		all_med_devs.push_back(fabs(alignment_median_dev - all_abs_devs[i]));
	}
	return MAD_CONSTANT * (median(all_med_devs));
}

double SQLiteConstructor::calculate_MAD_quicktree_sample(vector<Sequence> * inseqs, vector<Sequence> * inuserseqs) {
	srand(time(NULL));
	vector<int> rands;
	vector<Sequence> tseqs;
	vector<Sequence> tseqs2;
	if (inseqs->size() > 0) {
		for (int i = 0; i < RANDNUM; i++) {
			int n = rand() % inseqs->size();
			bool x = false;
			for (int j = 0; j < rands.size(); j++) {
				if (n == rands[j]) {
					x = true;
				}
				continue;
			}
			if (x == true) {
				i--;
			} else {
				rands.push_back(n);
			}
		}
		sort(rands.begin(), rands.end());
		for (int i = 0; i < RANDNUM; i++) {
			tseqs.push_back(inseqs->at(rands[i]));
		}
	}
	if (inuserseqs->size() > 0) {
		rands.clear();
		for (int i = 0; i < RANDNUM; i++) {
			int n = rand() % inuserseqs->size();
			bool x = false;
			for (int j = 0; j < rands.size(); j++) {
				if (n == rands[j]) {
					x = true;
				}
				continue;
			}
			if (x == true) {
				i--;
			} else {
				rands.push_back(n);
			}
		}
		sort(rands.begin(), rands.end());
		for (int i = 0; i < RANDNUM; i++) {
			tseqs2.push_back(inuserseqs->at(rands[i]));
		}
	}
	make_mafft_multiple_alignment(&tseqs, &tseqs2);
	return calculate_MAD_quicktree();
}

void SQLiteConstructor::load_info_on_immediate_children_of_taxon_id_into(string parent_taxon_id, vector<string> & child_ncbi_tax_ids,
		vector<string> & child_tax_names) {

	/* just get the immediate children (depth 1) from the parent taxon, and
	 * load their ncbi ids and taxon names into the passed vectors
	 */

	Database conn(db);
	Query query(conn);
	string sql = "SELECT ncbi_id, name FROM taxonomy WHERE parent_ncbi_id = ";
	sql += parent_taxon_id;
	sql += " and name_class = 'scientific name' and name_class not like '%environmental%' and name not like '%environmental%'";
	query.get_result(sql);

	while (query.fetch_row()) {
		string this_child_id = query.getstr();
		string this_child_name = query.getstr();
		child_ncbi_tax_ids.push_back(this_child_id);
		child_tax_names.push_back(this_child_name);
	}
	query.free_result();
}

void SQLiteConstructor::saturation_tests(vector<string> taxon_ids_to_be_tested, vector<string> taxon_names_to_be_tested,
		vector<Sequence> * source_db_seqs_to_assign) {

	/* The saturation tests use the mad score to subdivide the set of all taxa into the smaller sets, which will
	 * be used to make the original alignments that we will use to begin the profiling process.
	 *
	 * IF we are using the ncbi taxonomy as the guide tree for the tests:
	 *
	 * We start with the ncbi taxon ids in the passed taxon_ids_to_be_tested vector, which at this point, as far
	 * as i can tell, is just the search clade from the config file. for every one of these taxon ids, we pull
	 * all the sequences from keep_seqs that represent children of this taxon, using the
	 * filter_seqs_against_valid_children_of_ncbi_taxon_id function, align these sequences with mafft, and feed
	 * this alignment to quicktree to create a distance matrix. from this distance matrix, we calculate the mad
	 * score for the taxon. If the mad score passes, then we record the set of all sequences for the current taxon into an alignment.
	 *
	 * Things are different if we are using a user-supplied guide tree... or if we are doing an update...
	 *
	 */

	/* sas:
	 * changed this to accept a set of taxon_ids_to_be_tested and names
	 * this should allow for more flexible input when updating
	 * species
	 *
	 * the standard run will simply put one name and one id
	 * in the vectors
	 */

	vector<Sequence> seqs_to_be_assigned;

	// add seqs from the source db to the working set
	int n_source_db_seqs = source_db_seqs_to_assign->size();
	for (int i = 0; i < n_source_db_seqs; i++) {
		seqs_to_be_assigned.push_back(source_db_seqs_to_assign->at(i));
	}

	// add user seqs to the working set
	if (userfasta == true) {
		int n_user_seqs = user_seqs->size();
		for (int i = 0; i < n_user_seqs; i++) {
			seqs_to_be_assigned.push_back(user_seqs->at(i));
		}
	}

	// remember the alignments that we add
	vector<int> original_alignments_added;
	vector<string> seq_set_filenames; // here they will be not node names but numbers <- what numbers? this is unclear

	if (ncbi_saturation == true) {
		string tax_name_to_test;
		string tax_id_to_test;

		while (!taxon_names_to_be_tested.empty()) {

			/* each iteration of the saturation test takes the next taxon off the stack
			 * and checks (using the mad score) if it needs to be subdivided.
			 *
			 * If it does then we check to first check to make sure we're not subdividing
			 * a genus into a multitude of single-species alignments. If we are (happens
			 * fairly frequently with large genera and causes dramatic slowdowns) then we
			 * just keep the genus as a unit anyway and add it to the alignments. If it
			 * isn't a genus getting broken into lots of species though, then we subdivide,
			 * add the taxon's children to the stack, and continue. We keep running the
			 * tests until all taxa (or taxa they have been subdivided into) have passed
			 * the MAD test.
			 *
			 * sequences in seqs_to_align that haven't been assigned at the end of the
			 * tests are treated by the code that follows, that deals with leftovers.
			 */

			// first pop this search taxon off the stack
			tax_id_to_test = taxon_ids_to_be_tested.back();
			taxon_ids_to_be_tested.pop_back();
			tax_name_to_test = taxon_names_to_be_tested.back();
			taxon_names_to_be_tested.pop_back();

			// initialize containers to hold the child sequences of this taxon
//			vector<Sequence> * test_tax_db_child_seqs = new vector<Sequence>();
//			vector<Sequence> * test_tax_user_child_seqs = new vector<Sequence>();
//			test_tax_db_child_seqs->empty();
//			test_tax_user_child_seqs->empty();

			vector<Sequence> * test_tax_child_seqs = new vector<Sequence>();
			test_tax_child_seqs->empty();

			// get all the child seqs
			get_child_seqs_for_ncbi_tax_id(tax_id_to_test, source_db_seqs_to_assign, test_tax_child_seqs);
//			get_seqs_for_names_user(tax_id_to_test, test_tax_user_child_seqs); // TODO: still need better function names

			// precalculate counts of seqs to process
//			int n_db_seqs_to_test = test_tax_db_child_seqs->size();
//			int n_user_seqs_to_test = test_tax_user_child_seqs->size();
//			int n_seqs_to_test = test_tax_child_seqs->size();

			if (test_tax_child_seqs->size() == 0) {
				continue; // nothing to test for this taxon
			}

//			cout << "\ntaxon " << tax_name_to_test << " (ncbi id = " << tax_id_to_test << "): " << test_tax_child_seqs->size() << " seqs" << endl; //<< test_tax_user_child_seqs->size() << " user seqs; unassigned seqs remaining: " << seqs_to_be_assigned.size() << endl;
			if (test_tax_child_seqs->size() == 1) {

				// there can be no gaps in a one-sequence alignment, so just find the seq and set it as aligned
//				if (n_db_seqs_to_test > 0) {

				test_tax_child_seqs->at(0).set_aligned_sequence(test_tax_child_seqs->at(0).get_sequence());
				// TODO: should use the ncbi taxon id or the gi to find/remove db sequences
				remove_all_matching_seqs_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_child_seqs->at(0).get_taxon_name());

//				} else {
//					test_tax_user_child_seqs->at(0).set_aligned_sequence(test_tax_user_child_seqs->at(0).get_sequence());
//					remove_seq_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_user_child_seqs->at(0).get_taxon_name());
//				}

				// add the alignment record to the db
//				int alignid = gene_db.add_original_alignment(tax_id_to_test, test_tax_db_child_seqs, test_tax_user_child_seqs);
				int alignid = gene_db.add_original_alignment(tax_id_to_test, test_tax_child_seqs);
				if (updateDB == true)
					gene_db.toggle_alignment_update(alignid);
				original_alignments_added.push_back(alignid);

			} else { // there are multiple sequences; need to make an alignment
				double mad;
				cout << "using mafft for preliminary alignment" << endl;

				if (test_tax_child_seqs->size() < 2) {
					// can't calculate mad for fewer than three seqs; just do alignment
//					make_mafft_multiple_alignment(test_tax_db_child_seqs, test_tax_user_child_seqs);
					make_mafft_multiple_alignment(test_tax_child_seqs);
					mad = 0;

				} else if (test_tax_child_seqs->size() < MAX_SEQS_PER_ALIGNMENT) {
					// as long as we don't have too many seqs, do the alignment and get the mad score
//					make_mafft_multiple_alignment(test_tax_db_child_seqs, test_tax_user_child_seqs);
					make_mafft_multiple_alignment(test_tax_child_seqs);
					mad = calculate_MAD_quicktree();

				} else {
					// if there are too many seqs, make sure this alignment gets broken up
					mad = mad_cutoff + 1;
				}

				cout << "mad: " << mad << endl; // TODO: fix this so that if the alignment is too large it says that

				// check the mad score and decide what to do with the alignment
				bool add_alignment = false;

				// if the mad score is good
				if (mad <= mad_cutoff) {
					add_alignment = true;

				// if the mad score is bad
				} else {

					// check to make sure we're not about to dismantle a genus into many 1-species alns
					bool explode_taxon;

					Database conn(db);
					Query query(conn);

					// get the rank of the current taxon
					string sql = "select node_rank from taxonomy where name_class == 'scientific name' and ncbi_id == " + tax_id_to_test;
					query.get_result(sql);

					string rank = "";
					while (query.fetch_row())
						rank = query.getstr();

					// if it's not a genus, blow it up
					if (rank != "genus") {
						explode_taxon = true;

					} else { // if it is a genus...

						// get all the immediate taxonomic children of the genus
						vector<string> child_ids;
						vector<string> child_names; // just an argument required by the function; sloppy since we don't use this info
						load_info_on_immediate_children_of_taxon_id_into(tax_id_to_test, child_ids, child_names);

						// check how many are single species
						int n_single_taxon_children = 0;
						int n_children = 0;
						for (int i = 0; i < child_ids.size(); i++) {

							// first get lval/rval for this child
							string sql = "select left_value, right_value from taxonomy where ncbi_id == " + tax_id_to_test;
							query.get_result(sql);
							int lval, rval;
							while (query.fetch_row()) {
								lval = query.getval();
								rval = query.getval();
							}

							// increment counters
							n_children++;
							if ((rval - lval) < 2)
								// no children
								n_single_taxon_children++;
						}

						// only explode the genus if there is at least one multi-taxon child, and not too many singletons
						if (n_children > n_single_taxon_children && n_single_taxon_children <= MAX_ONE_TAXON_CHILD_ALNS) {
							explode_taxon = true;
						} else {
							explode_taxon = false;
							cout << "not splitting this genus as it would add " << n_single_taxon_children << " single-species alignments" << endl;
						}
					}

					if (explode_taxon)
						// push the children into taxon_ids_to_be_tested
						load_info_on_immediate_children_of_taxon_id_into(tax_id_to_test, taxon_ids_to_be_tested, taxon_names_to_be_tested);
					else
						add_alignment = true;
				}

				// if the mad score was good or we refused to break up a genus
				if (add_alignment) {

					// read in the seqs from the last alignment
//					update_seqs_using_last_alignment(test_tax_db_child_seqs, test_tax_user_child_seqs);
					update_seqs_using_last_alignment(test_tax_child_seqs);

					// remove the seqs we are assigning from the to-be-assigned vector
/*					for (int i = 0; i < test_tax_db_child_seqs->size(); i++) {
						remove_seq_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_db_child_seqs->at(i).get_taxon_name());
					}
					for (int i = 0; i < test_tax_user_child_seqs->size(); i++) {
//                        cout << "extracting user sequence from alignment by name: " << test_tax_user_child_seqs->at(i).get_taxon_name() << endl;
						remove_seq_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_user_child_seqs->at(i).get_taxon_name());
					} */


//					cout << endl << "removing sequences (currently there are " << seqs_to_be_assigned.size() << ")" << endl;
					for (int i = 0; i < test_tax_child_seqs->size(); i++) {
//						cout << test_tax_child_seqs->at(pp).get_taxon_name() << endl;
						remove_all_matching_seqs_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_child_seqs->at(i).get_taxon_name());
					}
//					cout << endl << "removed " << pp << "sequences, now there are " << seqs_to_be_assigned.size() << endl;

//					int alignid = gene_db.add_original_alignment(tax_id_to_test, test_tax_db_child_seqs, test_tax_user_child_seqs);
					int alignid = gene_db.add_original_alignment(tax_id_to_test, test_tax_child_seqs);
					if (updateDB == true) {
						gene_db.toggle_alignment_update(alignid);
                    }

					original_alignments_added.push_back(alignid);
				}
			}
			delete (test_tax_child_seqs);
			
			/*
			cout << endl << "remaining in seqs_to_be_assigned:" << endl;
			for (int q = 0; q < seqs_to_be_assigned.size(); q++) {
				cout << seqs_to_be_assigned.at(q).get_taxon_name() << endl;
			}
			cout << endl;
			*/

//			delete (test_tax_db_child_seqs);
//			delete (test_tax_user_child_seqs);
            
		} // END NCBI SATURATION
        
	} else { // user guide tree
		/*
		 * The idea here is to use the tree structure as the guide for the alignment and the
		 * breaking up of the groups. So the steps are
		 *
		 * 1) the stack is nodes
		 * 2) pop a node off and get all the seqs that are in that
		 * 3) do everything, if there is saturation push the children in there
		 * 3a) as part of the do everything bit, the mafft alignments should be recieving the node
		 *     as a guide tree for the alignment
		 */
		//for now ignoring the names that are sent because it is just the root, for update it should be the nodes
		cout << "using user guide tree" << endl;
		vector<Node *> nodes;
		if (updateDB == true) {
			for (int i = 0; i < taxon_names_to_be_tested.size(); i++) {
				bool found = false;
				for (int j = 0; j < userguidetree->getNodeCount(); j++) {
					if (userguidetree->getNode(j)->getName() == taxon_names_to_be_tested[i]) {
						nodes.push_back(userguidetree->getNode(j));
						found = true;
						break;
					}
				}
				if (found == false) {
					cerr << "problem updating and finding this node : " << taxon_names_to_be_tested[i] << endl;
					exit(0);
				}
			}
		} else {
			nodes.push_back(userguidetree->getRoot()); //this is different for update
		}
		while (!nodes.empty()) {
			Node * curnode = nodes.back();
			nodes.pop_back();
			vector<Sequence> * test_tax_db_seqs = new vector<Sequence>();
			vector<Sequence> * test_tax_user_child_seqs = new vector<Sequence>();
			test_tax_user_child_seqs->empty();
			test_tax_db_seqs->empty();
			get_seqs_for_nodes(curnode, source_db_seqs_to_assign, test_tax_db_seqs);
			get_seqs_for_user_nodes(curnode, test_tax_user_child_seqs);
			cout << test_tax_db_seqs->size() << " " << test_tax_user_child_seqs->size() << endl;
			if (test_tax_db_seqs->size() + test_tax_user_child_seqs->size() == 1) { //just one sequence in the group
				cout << curnode->getName() << " " << test_tax_db_seqs->size() << " " << test_tax_user_child_seqs->size() << endl;
				//make file
				for (int i = 0; i < test_tax_db_seqs->size(); i++) {
					test_tax_db_seqs->at(i).set_aligned_sequence(test_tax_db_seqs->at(i).get_sequence());
					remove_all_matching_seqs_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_db_seqs->at(i).get_taxon_name());
				}
				//user seqs
				for (int i = 0; i < test_tax_user_child_seqs->size(); i++) {
					test_tax_user_child_seqs->at(i).set_aligned_sequence(test_tax_user_child_seqs->at(i).get_sequence());
					remove_all_matching_seqs_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_user_child_seqs->at(i).get_taxon_name());
				}
				string alignment_name = curnode->getName();
				
				

				// NOTE: this is getting deactivated for now, because user seqs are no longer treated separately and
				// the add_original_alignment function only accepts one argument. THIS WILL BREAK THINGS IF TRYING TO RUN THE UPDATE.
//				int alignid = gene_db.add_original_alignment(alignment_name, test_tax_db_seqs, test_tax_user_child_seqs);
				int alignid = 0; // FAKE. JUST HERE AS A PLACEHOLDER



				if (updateDB == true)
					gene_db.toggle_alignment_update(alignid);
				original_alignments_added.push_back(alignid);
			} else if (test_tax_db_seqs->size() + test_tax_user_child_seqs->size() == 0) {
				continue;
			} else {
				/*
				 * multiple sequences
				 */
				cout << curnode->getName() << " " << test_tax_db_seqs->size() << endl;
				double mad;
				if (test_tax_db_seqs->size() + test_tax_user_child_seqs->size() > 2) {
					if (test_tax_db_seqs->size() + test_tax_user_child_seqs->size() < MAX_SEQS_PER_ALIGNMENT) {
						// TODO: add input tree for mafft
						make_mafft_multiple_alignment(test_tax_db_seqs, test_tax_user_child_seqs);
						mad = calculate_MAD_quicktree();
						/*}else if(test_tax_db_seqs->size() +temp_user_seqs->size()< 10000){
						 //need to make this happen 10 tens and average
						 mad = 0;
						 for (int i=0;i<10;i++)
						 mad = mad + (calculate_MAD_quicktree_sample(test_tax_db_seqs,temp_user_seqs)/10.0);
						 mad = mad * 2; //make sure is conservative*/
					} else { //if it is really big
						mad = mad_cutoff + 1; //make sure it gets broken up
					}
				} else {
					make_mafft_multiple_alignment(test_tax_db_seqs, test_tax_user_child_seqs);
					mad = 0;
				}
				cout << "mad: " << mad << endl;
				//if mad scores are good, store result
				if (mad <= mad_cutoff) {
					
					
					
					// NOTE: this is getting deactivated for now, because user seqs are no longer treated separately and
					// the add_original_alignment function only accepts one argument. THIS WILL BREAK THINGS IF TRYING TO RUN THE UPDATE.
//					update_seqs_using_last_alignment(test_tax_db_seqs, test_tax_user_child_seqs);
					


					for (int i = 0; i < test_tax_db_seqs->size(); i++) {
						remove_all_matching_seqs_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_db_seqs->at(i).get_taxon_name());
					}
					for (int i = 0; i < test_tax_user_child_seqs->size(); i++) {
						remove_all_matching_seqs_from_vector_by_taxon_name(&seqs_to_be_assigned, test_tax_user_child_seqs->at(i).get_taxon_name());
					}
					string alignment_name = curnode->getName();
					
					
					
					
					// NOTE: this is getting deactivated for now, because user seqs are no longer treated separately and
					// the add_original_alignment function only accepts one argument. THIS WILL BREAK THINGS IF TRYING TO RUN THE UPDATE.
//					int alignid = gene_db.add_original_alignment(alignment_name, test_tax_db_seqs, test_tax_user_child_seqs);
					int alignid = 0; // FAKE. JUST HERE AS A PLACEHOLDER




					if (updateDB == true)
						gene_db.toggle_alignment_update(alignid);
					original_alignments_added.push_back(alignid);
				}
				//if mad scores are bad push the children into names
				else {
					//need to get the children
					for (int i = 0; i < curnode->getChildCount(); i++) {
						nodes.push_back(curnode->getChild(i));
					}
				}
			}
			delete (test_tax_db_seqs);
			delete (test_tax_user_child_seqs);
		}
	}
	/*
	 * deal with the singletons
	 *
	 * singletons should be sequences that either don't have any data in the tree
	 * that is input or in the ncbi database if that is the tree to be used
	 */

	int n_leftovers = seqs_to_be_assigned.size();
	cout << "\n" << n_leftovers << " unassigned seqs";

	if (assignleftovers) {

		cout << ". picking where these should go (somewhat experimental)" << endl;

		// open the source db to get info about the leftovers
		SQLiteDBController dbc = SQLiteDBController(db);

		for (int i = 0; i < n_leftovers; i++) {

			// make a sequence object for this leftover seq
			Sequence this_leftover = seqs_to_be_assigned[i];

			// find the alignment for which this leftover seq has the highest affinity
			int best_score = 0;
			int best_match = 0;
			for (int j = 0; j < original_alignments_added.size(); j++) {
				vector<Sequence> tempseqs;
				bool aligned = true;
				gene_db.load_unaligned_seqs_from_original_alignment_into(tempseqs, original_alignments_added[j]);
				int tscore = get_single_to_group_seq_score(this_leftover, tempseqs);
				if (tscore > best_score) {
					best_score = tscore;
					best_match = j;
				}
			}

			// get info about the best matching alignment for this leftover
			int best_aln_for_this_leftover = original_alignments_added[best_match];
			string best_aln_db_name = gene_db.get_original_alignment_name_for_db_id(best_aln_for_this_leftover);
			string best_aln_tax_name = dbc.get_sci_name_for_ncbi_tax_id(atoi(best_aln_db_name.c_str()));

			cout << "adding leftover seq for " << this_leftover.get_taxon_name();
			if (this_leftover.is_user_seq()) {
				cout << " (user sequence)";
			} else {
				cout << " (gi " << this_leftover.get_ncbi_gi_number() << ")";
            }
			cout << " to " + best_aln_tax_name << endl;

			// get all the seqs for the best matching alignment and make a new alignment including this leftover
			vector<Sequence> tempseqs;
			gene_db.load_unaligned_seqs_from_original_alignment_into(tempseqs, best_aln_for_this_leftover);
			tempseqs.push_back(this_leftover);
			make_mafft_multiple_alignment(&tempseqs);

			// get the newly aligned seq, and add this aligned seq to the alignment in the db
			retrieve_aligned_sequence_from_last_alignment_for_seq(&this_leftover);
			gene_db.associate_sequence_with_alignment(best_aln_for_this_leftover, this_leftover);

			// read mafft
			vector<Sequence> aligned_seqs;
			load_sequences_from_last_alignment_into(aligned_seqs);

			gene_db.update_align_seqs(best_aln_for_this_leftover, aligned_seqs);
            
/*            if (i > 5) {
                exit(0);
            } */

			if (updateDB == true)
				gene_db.toggle_alignment_update(original_alignments_added[best_match]);
		}

	} else { // do not assign leftovers
		if (n_leftovers > 0) {
			cout
					<< " (written to the '.leftovers' file) will not be included in the run.\n(to include them, issue the keyword 'assignleftovers' in the config file, but be aware \nthey may be assigned in surprising places!)."
					<< endl;

			// write the leftovers to a file
			ofstream leftfile;
			leftfile.open((gene_name + ".leftovers").c_str(), ios::out);
			leftfile << "ncbi_tax_id\ttaxon_name\tgi\tsource" << endl;
			for (int i = 0; i < n_leftovers; i++) {
				Sequence this_leftover = seqs_to_be_assigned[i];
				leftfile << this_leftover.get_ncbi_tax_id() << "\t" << this_leftover.get_taxon_name() << "\t" << this_leftover.get_ncbi_gi_number() << "\t"
						<< this_leftover.get_source() << endl;
			}
			leftfile.close();
		}
	}
	cout << "\nfinished assembly" << endl;
}

/*
 * comparing a sequences to a group of sequences and returning the best score
 * 
 * this can be helpful when deciding with which group a sequence goes
 *
 * ASSUMPTIONS: the sequences are all pointing the right direction
 */
int SQLiteConstructor::get_single_to_group_seq_score(Sequence & inseq, vector<Sequence> & ginseqs) {
	vector<int> scores;
	SBMatrix mat = swps3_readSBMatrix("EDNAFULL");
	for (int i = 0; i < ginseqs.size(); i++) {
		double maxide = 0;
		int ret = swps3_maxscores(mat, &inseq, &ginseqs[i]);
		double tsc = double(ret);
		//cout <<i << " " << j << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << endl;
		if (std::numeric_limits<double>::infinity() != tsc) {
			scores.push_back(tsc);
		}
	}
	int maxide = 0;
	for (int i = 0; i < scores.size(); i++) {
		if (scores[i] > maxide)
			maxide = scores[i];
	}
}

/*
 * this stores the gi numbers for reference
 */

void SQLiteConstructor::write_gi_numbers(vector<Sequence> * dbs) {
	for (int i = 0; i < dbs->size(); i++) {
		//gifile << dbs->at(i).get_tax_id() << "\t"; //don't need this anymore
		gifile << dbs->at(i).get_ncbi_tax_id() << "\t";
		//gifile << dbs->at(i).get_accession() << endl;
		gifile << dbs->at(i).get_ncbi_gi_number() << "\t";
		gifile << dbs->at(i).get_taxon_name() << endl;
	}
}

/*
 * this stores the names for the user seqs, mostly this is for updates
 */
void SQLiteConstructor::write_user_numbers() {
	for (int i = 0; i < user_seqs->size(); i++) {
		ufafile << user_seqs->at(i).get_taxon_name() << "\t";
		ufafile << user_seqs->at(i).get_ncbi_tax_id() << "\t";
		ufafile << endl;
	}
}

/*
 * this is primarily used for add the seqs from a file to the dbseq and rc vectors 
 * most for updating alignments
 *
 * all the sequences need to be in the database for this to work
 */
void SQLiteConstructor::add_seqs_from_db_to_seqs_vector(string alignment_name, vector<Sequence> * keep_seqs, vector<Sequence> & storedseqs) {
	FastaUtil fu;
	vector<Sequence> tseqs;
	int alignment_id = gene_db.get_original_alignment_id_by_name(alignment_name);
	gene_db.load_unaligned_seqs_from_original_alignment_into(tseqs, alignment_id);
	cout << "seqs from " << alignment_name << ": " << tseqs.size() << endl;
	Database conn(db);
	for (unsigned int i = 0; i < tseqs.size(); i++) {

		/*		size_t found;
		 //if it is a user seq
		 found = tseqs[i].get_label().find("user_"); */

		if (tseqs[i].is_user_seq()) {

			//TODO: not sure if I need this
			if (usertree == true) {
				cout << "matching user fasta seqs to user guide tree" << endl;
				int count = 0;
				for (int i = 0; i < userguidetree->getExternalNodeCount(); i++) {
					string tname = userguidetree->getExternalNode(i)->getName();
					cout << tname << endl;
					if (tname == tseqs[i].get_taxon_name() || tname == tseqs[i].get_ncbi_tax_id() || ("user_" + tname == tseqs[i].get_taxon_name())
							|| ("user_" + tname == tseqs[i].get_ncbi_tax_id())) {
						user_fasta_node_map[&tseqs[i]] = userguidetree->getExternalNode(i);
						count += 1;
					}
				}
				cout << "matches: " << count << " prop:" << count / user_seqs->size() << endl;
			}
			keep_seqs->push_back(tseqs[i]);
		} else { //it is a dbseq
			keep_seqs->push_back(tseqs[i]);
		}
	}
}

double SQLiteConstructor::get_usertree_keepseq_overlap(vector<Sequence> * keep_seqs) {
	set<string> tree_names;
	for (int i = 0; i < userguidetree->getExternalNodeCount(); i++) {
		tree_names.insert(userguidetree->getExternalNode(i)->getComment());
	}
	double ccount = 0;
	for (int i = 0; i < keep_seqs->size(); i++) {
		//need to change this to be hierachical
		if (tree_names.count(keep_seqs->at(i).get_ncbi_tax_id()) == 1)
			ccount += 1;
	}
	return ccount / (double) keep_seqs->size();
}

Tree * SQLiteConstructor::get_user_guide_tree_obj() {
	return userguidetree;
}

bool SQLiteConstructor::get_updatestatus() {
	return updateDB;
}

string SQLiteConstructor::get_genedb() {
	return gene_db_name;
}

//void SQLiteConstructor::remove_seq_from_vector_by_taxon_name(vector<Sequence> * v, string taxon_name) {
void SQLiteConstructor::remove_all_matching_seqs_from_vector_by_taxon_name(vector<Sequence> * v, string taxon_name) {
//	int position = 0;
	vector<int> matched_positions;
//	int sz = v->size();
	for (int i = 0; i < v->size(); i++) {
		if (taxon_name == (*v)[i].get_taxon_name()) {
			// decrement the index by the number of previous matches because they will be deleted before we get to this one
			matched_positions.push_back(i - matched_positions.size());
//			break;
		}
	}

	for (int i = 0; i < matched_positions.size(); i++) {
		v->erase(v->begin() + matched_positions.at(i));
	}

//	if (v->size() != sz - 1) {
	if (matched_positions.size() < 1) {
		cerr << "sequence '" + taxon_name + "' could not be removed from the vector. it may not have been found" << endl;
		exit(1);
//	} else {
//		cout << "erased " << matched_positions.size() << " matching sequences" << endl;
	}
}

//void SQLiteConstructor::remove_seq_from_vector_by_ncbi_id(vector<Sequence> * v, string ncbi_id) {
void SQLiteConstructor::remove_all_matching_seqs_from_vector_by_ncbi_id(vector<Sequence> * v, string ncbi_id) {
//	int position = 0;
	vector<int> matched_positions;
//	int sz = v->size();
	for (int i = 0; i < v->size(); i++) {
		if (ncbi_id == (*v)[i].get_ncbi_tax_id()) {
			// decrement the index by the number of previous matches because they will be deleted before we get to this one
			matched_positions.push_back(i - matched_positions.size());
//			position = i;
//			break;
		}
	}
//	v->erase(v->begin() + position);
	for (int i = 0; i < matched_positions.size(); i++) {
		v->erase(v->begin() + matched_positions.at(i));
	}

	if (matched_positions.size() < 1) {
		cerr << "sequence '" + ncbi_id + "' could not be removed from the vector. it may not have been found" << endl;
		exit(1);
	} else {
		cout << "erased " << matched_positions.size() << " matching sequences" << endl;
	}

//	if (v->size() != sz - 1) {
//		cerr << "sequence for ncbi taxon '" + ncbi_id + "' could not be removed from the vector. it may not have been found" << endl;
//		exit(1);
//	}
}

//void SQLiteConstructor::update_seqs_using_last_alignment(vector<Sequence> * db_seqs_to_update, vector<Sequence> * user_seqs_to_update) {
void SQLiteConstructor::update_seqs_using_last_alignment(vector<Sequence> * seqs_to_update) {
	FastaUtil fu;
	vector<Sequence> aligned_seqs;
	string outfile = genefoldername + "outfile";
	fu.read_aligned_fasta_into(aligned_seqs, outfile); // look here next...
	for (int i = 0; i < aligned_seqs.size(); i++) {

		bool set = false;

		for (int j = 0; j < seqs_to_update->size(); j++) {
			if (aligned_seqs[i].get_ncbi_tax_id() == seqs_to_update->at(j).get_ncbi_tax_id()) {

//                	cout << "updating db sequence: " << aligned_seqs[i].get_ncbi_tax_id() << endl;

				seqs_to_update->at(j).set_aligned_sequence(aligned_seqs[i].get_sequence());
//                    cout << "added to db seqs" << endl;
				set = true;
				break;
			}
		}
        
//        cout << "adding to incoming sequences: " << aligned_seqs.at(i).get_taxon_name() << "... ";
        
        /*
		for (int j = 0; j < db_seqs_to_update->size(); j++) {
			if (aligned_seqs[i].get_ncbi_tax_id() == db_seqs_to_update->at(j).get_ncbi_tax_id()) {
				db_seqs_to_update->at(j).set_aligned_sequence(aligned_seqs[i].get_sequence());
                cout << "added to db seqs" << endl;
				set = true;
				break;
			}
		}
		if (set == false) {
			for (int j = 0; j < user_seqs_to_update->size(); j++) {
				if (aligned_seqs[i].get_taxon_name() == user_seqs_to_update->at(j).get_taxon_name()) {
					user_seqs_to_update->at(j).set_aligned_sequence(aligned_seqs[i].get_sequence());
                    cout << "added to user seqs" << endl;
					set = true;
					break;
				}
			}
		} */

/*
        if (aligned_seqs[i].is_user_seq() == true) {
            for (int j = 0; j < user_seqs_to_update->size(); j++) {
				if (aligned_seqs[i].get_taxon_name() == user_seqs_to_update->at(j).get_taxon_name()) {

//		        	cout << "updating user sequence: " << aligned_seqs[i].get_taxon_name() << endl;

					user_seqs_to_update->at(j).set_aligned_sequence(aligned_seqs[i].get_sequence());
//                    cout << "added to user seqs" << endl;
                    set = true;
					break;
				}
			}
        } else { // not a user sequence
            for (int j = 0; j < db_seqs_to_update->size(); j++) {
                if (aligned_seqs[i].get_ncbi_tax_id() == db_seqs_to_update->at(j).get_ncbi_tax_id()) {

//                	cout << "updating db sequence: " << aligned_seqs[i].get_ncbi_tax_id() << endl;

                    db_seqs_to_update->at(j).set_aligned_sequence(aligned_seqs[i].get_sequence());
//                    cout << "added to db seqs" << endl;
                    set = true;
                    break;
                }
            }
        }
        */
        
		if (set == false) {
			cout << "error, could not find match for " << aligned_seqs[i].get_taxon_name() << ". phlawd will exit" << endl;
			exit(1);
		}
	}
}

void SQLiteConstructor::retrieve_aligned_sequence_from_last_alignment_for_seq(Sequence * seq) {

	/* checks the last mafft alignment to see if the passed sequence is in the
	 * alignment (looks for the sequence's label), and sets the aligned sequence
	 * for the passed sequence object to the sequence from the alignment.
	 */

	// read in the last alignment
	FastaUtil fu;
	vector<Sequence> tempalseqs;
	string outfile = genefoldername + "outfile";
	fu.read_aligned_fasta_into(tempalseqs, outfile);

	// look for the seq
	bool seq_was_found = false;
	for (int i = 0; i < tempalseqs.size(); i++) {

		// if we find it, get the aligned seq from the file and update the seq
		if (tempalseqs[i].is_user_seq()) {
			if (tempalseqs[i].get_taxon_name() == seq->get_taxon_name()) {
				seq->set_aligned_sequence(tempalseqs[i].get_sequence());
				seq_was_found = true;
				break;
			}
		} else if (tempalseqs[i].get_ncbi_tax_id() == seq->get_ncbi_tax_id()) {
			seq->set_aligned_sequence(tempalseqs[i].get_sequence());
			seq_was_found = true;
			break;
		}
	}

	if (!seq_was_found) { // if it wasn't there, something is probably wrong
		cout << "error, sequence '" << seq->get_taxon_name() << "' was not in the last alignment. phlawd will exit" << endl;
		exit(0);
	}
}

void SQLiteConstructor::load_sequences_from_last_alignment_into(vector<Sequence> & seqs) {

	FastaUtil fu;
//	vector<Sequence> tempalseqs;
	string outfile = genefoldername + "outfile";
	fu.read_aligned_fasta_into(seqs, outfile);

	/*	this is just duplicting effort
	 for (int i = 0; i < tempalseqs.size(); i++) {
	 tempalseqs[i].set_label(tempalseqs[i].get_label());
	 tempalseqs[i].set_aligned_sequence(tempalseqs[i].get_sequence());
	 temp_seqs->push_back(tempalseqs[i]);
	 } */
}
