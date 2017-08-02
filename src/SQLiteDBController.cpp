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
 * SQLiteDBController.h
 */

// TODO: rework this class to provide all interactivity with the source db, so other classes don't have to.
#include <string>
#include <map>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <cstring>
#include <sstream>
#include <stdio.h>
#include <regex.h>
#include <unistd.h>

using namespace std;

#include "libsqlitewrapped.h"
#include <sqlite3.h>
#include "GenBankReader.h"
#include "SQLiteDBController.h"
#include "utils.h"

SQLiteDBController::SQLiteDBController(string dbn) :
		db_name(dbn), division(""),
				count(0) {
}

template<class T>
inline std::string to_string(const T& t) {
	std::stringstream ss;
	ss << t;
	return ss.str();
}

string SQLiteDBController::get_sci_name_for_ncbi_tax_id(int ncbi_tax_id) {

	Database conn(db_name);
	Query query(conn);

	string sql = "select name from taxonomy where name_class == 'scientific name' and ncbi_id == ";
	sql += to_string(ncbi_tax_id);
	query.get_result(sql);

	string taxon_name = "";

	while (query.fetch_row())
		taxon_name = query.getstr();

	return taxon_name;
}

string SQLiteDBController::get_rank_for_ncbi_tax_id(int ncbi_tax_id) {

	Database conn(db_name);
	Query query(conn);

	string sql = "select node_rank from taxonomy where name_class == 'scientific name' and ncbi_id == ";
	sql += to_string(ncbi_tax_id);
	query.get_result(sql);

	string rank = "";

	while (query.fetch_row())
		rank = query.getstr();

	return rank;
}

bool SQLiteDBController::initiate() {
	bool ret = true;
	//check to see if the database exists
	ifstream ifile(db_name.c_str());
	if (ifile) {
		cout << "the file: " + db_name + " seems to exist" << endl;
		return false;
	}
	Database conn(db_name);
	Query query(conn);
	query.get_result(
			"create table taxonomy (id INTEGER PRIMARY KEY,ncbi_id INTEGER,name TEXT,name_class TEXT,node_rank TEXT,parent_ncbi_id INTEGER,edited_name TEXT,left_value INTEGER,right_value INTEGER);");
	query.free_result();
	query.get_result("CREATE INDEX taxonomy_left_value on taxonomy(left_value);");
	query.free_result();
	query.get_result("CREATE INDEX taxonomy_name on taxonomy(name);");
	query.free_result();
	query.get_result("CREATE INDEX taxonomy_ncbi_id on taxonomy(ncbi_id);");
	query.free_result();
	query.get_result("CREATE INDEX taxonomy_parent_ncbi_id on taxonomy(parent_ncbi_id);");
	query.free_result();
	query.get_result("CREATE INDEX taxonomy_right_value on taxonomy(right_value);");
	query.free_result();
	query.get_result("CREATE INDEX taxonomy_edited_name on taxonomy(edited_name);");
	query.free_result();

	query.get_result(
			"create table sequence (id INTEGER PRIMARY KEY,ncbi_id INTEGER,accession_id TEXT,identifier TEXT,description TEXT,seq LONGTEXT);");
	query.free_result();
	query.get_result("CREATE INDEX sequence_ncbi_id on sequence(ncbi_id);");
	query.free_result();
	query.get_result("CREATE INDEX sequence_accession_id on sequence(accession_id);");
	query.free_result();
	query.get_result("CREATE INDEX sequence_identifier on sequence(identifier);");
	query.free_result();

	query.get_result("create table information (id INTEGER PRIMARY KEY, name TEXT, value TEXT);");
	query.free_result();
	return ret;
}

static int callback(void *NotUsed, int argc, char **argv, char **azColName) {
	int i;
	for (i = 0; i < argc; i++) {
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	return 0;
}

string SQLiteDBController::create_name(string & tfilen) {
	size_t found;
	found = tfilen.find("\"");
	while (found != string::npos) {
		tfilen.replace(found, 1, "'");
		found = tfilen.find("\"", found + 1);
	}
	return tfilen;
}

string SQLiteDBController::create_edited_name(string & tfilen) {
	size_t found;
	found = tfilen.find(" ");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find(" ", found + 2);
	}

    string badchars = "~`!@#$%^&*()+={[}]:;\"'<,>.?/|\\";
    
    for (int i = 0; i < badchars.length(); i++) {
        replaceChar(tfilen, badchars.at(i));
    } 

    /*
	//(take out the parenthetical stuff too)
	found = tfilen.find("(");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find("(", found + 1);
	}
	found = tfilen.find(")");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find(")", found + 1);
	}
	found = tfilen.find("\"");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find("\"", found + 1);
	}
	found = tfilen.find("'");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find("'", found + 1);
	}
	found = tfilen.find(".");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find(".", found + 1);
	}
	found = tfilen.find("&");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find("&", found + 1);
	}
	found = tfilen.find(",");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find(",", found + 1);
	}
	found = tfilen.find("\\");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find("\\", found + 1);
	}
	found = tfilen.find("/");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find("/", found + 1);
	}
	found = tfilen.find(";");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find(";", found + 1);
	}
    found = tfilen.find(":");
	while (found != string::npos) {
		tfilen.replace(found, 1, "_");
		found = tfilen.find(":", found + 1);
	} */
	return tfilen;
}

void SQLiteDBController::replaceChar(string & instring, char s) {
    size_t found = instring.find(s);
	while (found != string::npos) {
		instring.replace(found, 1, "_");
		found = instring.find(s, found + 1);
	}
}

/*ednm = spls[1].replace('\"',"_").strip()
 ednm = ednm.replace("\'","_")
 ednm = ednm.replace("\\","_")
 ednm = ednm.replace("/","_")
 ednm = ednm.replace("(","_")
 ednm = ednm.replace(")","_")
 ednm = ednm.replace(".","_")
 ednm = ednm.replace("&","_")
 ednm = ednm.replace(",","_")
 ednm = ednm.replace(" ","_")
 */
void SQLiteDBController::load_seqs(string div, string ref, bool downl) {
	cout << "loading taxonomy" << endl;
	if (downl == true) {
		const char * cmd = "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
		cout << "downloading with wget" << endl;
		system(cmd);
		cmd = "tar -xzvf taxdump.tar.gz";
		cout << "untaring" << endl;
		system(cmd);
	}
	//read the nodes.dmp
	map<string, string> rank;
	map<string, string> parent_id;
	ifstream infile("nodes.dmp", ios::in);
	string line;
	vector<string> tokens;
	while (getline(infile, line)) {
		string del("|");
		tokens.clear();
		Tokenize(line, tokens, del);
		if (tokens.size() > 1) {
			for (int i = 0; i < tokens.size(); i++) {
				TrimSpaces(tokens[i]);
			}
			string ncbi_id = tokens[0];
			rank[ncbi_id] = tokens[2];
			parent_id[ncbi_id] = tokens[1];
		}
	}
	infile.close();
	//read the names.dmp
	ifstream infile2("names.dmp", ios::in);
	sqlite3 *conn;
	int rc = sqlite3_open(db_name.c_str(), &conn);
	char *zErrMsg = 0;

	sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
	count = 0;
	while (getline(infile2, line)) {
		if (count % 100000 == 0) {
			cout << count << endl;
		}
		string del("|");
		tokens.clear();
		Tokenize(line, tokens, del);
		if (tokens.size() > 1) {
			for (int i = 0; i < tokens.size(); i++) {
				TrimSpaces(tokens[i]);
			}
			string gin = tokens[0];
			string nm = create_name(tokens[1]); //need to double quote the single quotes and backslash the quotes
			string nm_class = tokens[3];
			string ednm = create_edited_name(tokens[1]); //need to edit the names
			string sql = "insert into taxonomy (ncbi_id,name,name_class,node_rank,parent_ncbi_id,edited_name) values (";
			sql += gin + ",\"";
			sql += nm + "\",'";
			sql += nm_class + "','";
			sql += rank[gin] + "',";
			sql += parent_id[gin] + ",'";
			sql += ednm + "');";
			//query.execute(sql);
			rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
			//uncomment to get the names that don't commit, mostly bad quotes
//	    if (rc != 0)
//		cout << sql << endl;
		}
		count += 1;
	}
	sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
	infile2.close();
	sqlite3_close(conn);

	cout << "updating left/right values" << endl;
	Database cppconn(db_name);
	Query query(cppconn);
	string cmd = "select ncbi_id,parent_ncbi_id from taxonomy where name_class = 'scientific name';";
	query.get_result(cmd);
	while (query.fetch_row()) {
		int nc = query.getval();
		int pc = query.getval();
		if (parent_ncbi_map.count(pc) > 0) {
			parent_ncbi_map[pc].push_back(nc);
		} else {
			vector<int> tv;
			tv.push_back(nc);
			parent_ncbi_map[pc] = tv;

		}
	}
	//remove files
	remove("citations.dmp");
	remove("division.dmp");
	remove("gc.prt");
	remove("gencode.dmp");
	remove("readme.txt");

	//get the root and send to rebuild
	count = 0;
	rc = sqlite3_open(db_name.c_str(), &conn);
	sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
	rebuild_tree(1, 1, conn);
	sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
	sqlite3_close(conn);

	cout << "loading seqs" << endl;
	vector<string> filesToProcess;

	if (div.length() > 0) {
		division = div;
		vector<string> groups; // taxonomic groups to be downloaded
		if (division == "met" || division == "all") {
			groups.push_back("pri");
			groups.push_back("rod");
			groups.push_back("mam");
			groups.push_back("vrt");
			groups.push_back("inv");
			groups.push_back("bct");

			if (division == "all")
				groups.push_back("pln");

		} else {
			groups.push_back(division);
		}

		for (int i = 0; i < groups.size(); i++) {
			string cmd;
			string fnameString = "gb" + groups[i] + "*.seq.gz";
			if (downl == true) {
				cmd = "wget ftp://ftp.ncbi.nih.gov/genbank/" + fnameString ;
				cout << "downloading with wget" << endl;
				system(cmd.c_str());
				cmd = "gunzip -d gb" + groups[i] + "*.seq.gz";
				cout << "uncompressing" << endl;
				system(cmd.c_str());
//			}
			} else {
//			for (int i = 0; i < groups.size(); i++) {
//				cmd = "gunzip -d gb" + groups[i] + "*.seq.gz";
				cmd = "gunzip -d " + fnameString;
				cout << "uncompressing" << endl;
				system(cmd.c_str());
			}
		}

		// download daily updates
		if (downl == true) {
			string fnameString = "nc*.flat.gz";
			string cmd = "wget -nv ftp://ftp.ncbi.nih.gov/genbank/daily-nc/" + fnameString;
			cout << "downloading dailies with wget" << endl;
			system(cmd.c_str());
			cmd = "gunzip -d " + fnameString;
			cout << "uncompressing dailies" << endl;
			system(cmd.c_str());
		}

		// get the names of the files to use
		vector<string> file_names;
		cout << "getting file names for gb flat files" << endl;
		getdir(".", file_names);
		for (int i = 0; i < file_names.size(); i++) {
			for (int j = 0; j < groups.size(); j++) {
				if (file_names[i].find("gb" + groups[j]) != string::npos && file_names[i].substr(file_names[i].size() - 4, 4) == ".seq") {
					filesToProcess.push_back(file_names[i]);
				} else if (file_names[i].find("nc") != string::npos && file_names[i].substr(file_names[i].size() - 5, 5) == ".flat") {
					filesToProcess.push_back(file_names[i]);
				}
			}
		}
	}

	if (ref.length() > 0) { // if we're getting whole genomes

		refseq = ref;
		vector<string> groups; // taxonomic groups to be downloaded
		if (refseq == "metazoan" || refseq == "all") {
			groups.push_back("mitochondrion");
			groups.push_back("invertebrate");
			groups.push_back("vertebrate-mammalian");
			groups.push_back("vertebrate-other");

		} else if (refseq == "plant" || refseq == "all") {
			groups.push_back("plant");
			groups.push_back("plastid");

		} else if (refseq == "microbes" || refseq == "all") {
			groups.push_back("microbial");
			groups.push_back("plasmid");
			groups.push_back("protozoa");
		}

		if (refseq == "all") {
			groups.push_back("fungi");
			groups.push_back("viral");

		} else {
			groups.push_back(refseq);
		}

		for (int i = 0; i < groups.size(); i++) {
			string cmd;
			string fnameString = groups[i] + "*.genomic.gbff.gz";
			if (downl == true) {
				cout << "downloading with wget" << endl;
				cmd = "wget ftp://ftp.ncbi.nih.gov/refseq/release/" + groups[i] + "/" + fnameString; // ftp://ftp.ncbi.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
				system(cmd.c_str());
				cmd = "gunzip -d " + fnameString;
				cout << "uncompressing" << endl;
				system(cmd.c_str());
			//}
			} else {
//				for (int i = 0; i < groups.size(); i++) {
				cmd = "gunzip -d " + fnameString;
				cout << "uncompressing" << endl;
				system(cmd.c_str());
			}
		}

		// get the names of the files to use
		vector<string> file_names;
		cout << "getting file names for refseq flat files" << endl;
		getdir(".", file_names);
		for (int i = 0; i < file_names.size(); i++) {
			for (int j = 0; j < groups.size(); j++) {
				if (file_names[i].find(groups[j]) != string::npos && file_names[i].substr(file_names[i].size() - 5, 5) == ".gbff") {
					filesToProcess.push_back(file_names[i]);
//					cout << filen << endl;
//					GenBankReader gbr;
//					gbr.parse_file(filen, db_name);
//					remove(filen.c_str());
				}
			}
		}


	}

	// now process the files into the db
	for (int i = 0; i < filesToProcess.size(); i++) {
		cout << filesToProcess[i] << endl;
		GenBankReader gbr;
		gbr.parse_file(filesToProcess[i], db_name);
		remove(filesToProcess[i].c_str());
	}

	cout << "merging old names with new names" << endl;
	ifstream infile3("merged.dmp", ios::in);
	rc = sqlite3_open(db_name.c_str(), &conn);
	sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
	count = 0;
	while (getline(infile3, line)) {
		if (count % 100000 == 0) {
			cout << count << endl;
		}
		string del("|");
		tokens.clear();
		Tokenize(line, tokens, del);
		for (int i = 0; i < tokens.size(); i++) {
			TrimSpaces(tokens[i]);
		}
		if (tokens.size() > 1) {
			string sql = "update sequence set ncbi_id = ";
			sql += tokens[1];
			sql += " where ncbi_id = ";
			sql += tokens[0];
			sql += ";";
			rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
			if (rc != 0)
				cout << sql << endl;
		}
		count += 1;
	}
	infile3.close();
	sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
	sqlite3_close(conn);

	remove("merged.dmp");
	remove("names.dmp");
	remove("nodes.dmp");
	remove("delnodes.dmp");
	cout << "finished loading" << endl;
}

//private
int SQLiteDBController::rebuild_tree(int gid, int lft, sqlite3 * conn) {
	int rgt = lft + 1;
	vector<int> res;
	if (parent_ncbi_map.count(gid) > 0) {
		res = parent_ncbi_map[gid];
	}
	for (int i = 0; i < res.size(); i++) {
		if (res[i] == gid)
			continue;
		else
			rgt = rebuild_tree(res[i], rgt, conn);
	}
	string updcmd = "update taxonomy set left_value = " + to_string(lft) + ", right_value = " + to_string(rgt) + " where ncbi_id = " + to_string(gid) + ";"; //and name_class = scientific name
	//cout << updcmd << endl;
	//exit(0);
	int rc = sqlite3_exec(conn, updcmd.c_str(), 0, 0, 0);
	if (count % 100000 == 0)
		cout << count << endl;
	count += 1;
	return rgt + 1;
}
