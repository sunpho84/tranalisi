#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "node.hpp"

#include "fitter/fitter_parser.hpp"
#include "operations.hpp"

using namespace std;

void init_scanner();
void destroy_scanner();

// string s("a test string");

// typedef struct yy_buffer_state * YY_BUFFER_STATE;
// extern int yyparse();
// extern YY_BUFFER_STATE yy_scan_string(char * str);
// extern void yy_delete_buffer(YY_BUFFER_STATE buffer);

// node_t* parse(const string s)
// {
//   // yyscan_t scanner;
//   //   YY_BUFFER_STATE buf;
//   //   yylex_init(&scanner);
//   //   buf = yy_scan_string("replace me with the string youd like to scan", scanner);
//   //   yylex(scanner);
//   //   yy_delete_buffer(buf, scanner);
//   //   yylex_destroy(scanner);
// }

// int main()
// {
  
//   //   real val(10);
//   // uexp e(exp,&val);
//   // bexp s(sum,&val,&e);
//   // uexp ss(sin,&e);
//   // cout<<ss.eval()<<endl;
  
//   return 0;
// }

void fitter_lex_internal(char* buf,int &result,int max_size)
{
  result=0;
  
  // do
  //   {
  //     //point to last opened file
  //     ifstream& f=fin.back().first;
      
  //     //read
  //     f.get(*buf);
      
  //     //close file if end reached
  //     if(f.eof())
  // 	{
  // 	  print_input_files_stack();
  // 	  fin.pop_back();
  // 	}
  //     //check what read otherwise
  //     else
  // 	result=f.good();
  //   }
  // //exit loop when no more file can be read or result is valid
  // while(fin.size()>=1 and result==0);
}

int main()
{
  init_scanner();
  
  fitter_parser_parse(nullptr);
  
  destroy_scanner();
  
  return 0;
}

