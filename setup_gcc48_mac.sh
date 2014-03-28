# uncomment this line to install homebrew if necessary
#ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"

# use homebrew to install a gcc compiler with openmp support
brew tap homebrew/versions
brew install gcc48

# set up symlinks to gcc4-8 and g++4-8
ln -s /usr/local/bin/gcc-4.8 /usr/local/bin/gcc
ln -s /usr/local/bin/g++-4.8 /usr/local/bin/g++

# update PATH variables so new compilers are default 
cp ~/.bashrc bash_temp
printf "\nexport PATH=/usr/local/bin:\$PATH" >> bash_temp
printf "\nexport CPP=g++" >> bash_temp
printf "\nexport CXX=g++" >> bash_temp
printf "\nexport CC=gcc" >> bash_temp
mv bash_temp ~/.bashrc
source ~/.bashrc