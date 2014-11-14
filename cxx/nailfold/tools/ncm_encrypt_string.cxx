#include <vcl_iostream.h>

#include <nailfold/ncm_encrypt.h>

//
// Main function
int main( int argc, char* argv[] )
{
  vcl_ofstream ofs("u:/tmp/strings.txt");
  ofs << ncm_encrypt("your string here") << vcl_endl;
  ofs.close();

  return 0;
}


