#include <iostream>
#include <iostream>
#include<fstream>
#include "lexer_untyped.h"

using namespace std;

int main(int argc,  char** argv)
{
 //        auto filename=argv[1];
        auto filename="../macro_dr/test.txt";
        std::ifstream f(filename);
        std::string s;
        while (f)
          {
                std::string line;
                std::getline(f,line);
                s+=line+"\n";
          }
        std::cout<<"\ntest file \n"<<s<<"\n";
        auto p=dcli::extract_program(s);
        if(p){
                auto ss=(*p).str();
                std::cerr<<ss;
        }
        else
                std::cerr<<p.error()();
}
