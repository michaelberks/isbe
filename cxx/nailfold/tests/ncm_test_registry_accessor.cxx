#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_iomanip.h>
#include <vcl_ctime.h>

#include <nailfold/ncm_registry_accessor.h>

//
// Main function
int main( int argc, char* argv[] )
{
  ncm_registry_accessor ra;

  // Read value from an existing key:
  // HKEY_CURRENT_USER/Software/Microsoft/Notepad/lfFaceName (string)
  ra.get_current_user_key();
  ra.get_subkey("Software\\Microsoft\\Notepad");

  vcl_string s_value = "";
  ra.read_string("lfFaceName", s_value);


  // Create a new key
  ra.get_current_user_key();
  ra.create_subkey("Software\\University of Manchester\\Test Key");


  // Add a string name/value pair to new subkey
  ra.write_string("Name", "Phil Tresadern");

  // Read the name/value pair back
  ra.read_string("Name", s_value);
  vcl_cout << "Name = " << s_value << vcl_endl;


  // Add a numeric (DWORD) name/value pair to new subkey
  ra.write_numeric("Age", 32);

  // Read the name/value pair back
  int n_value = 0;
  ra.read_numeric("Age", n_value);
  vcl_cout << "Age = " << n_value << vcl_endl;


  // Delete created key
  ra.get_current_user_key();
  ra.get_subkey("Software\\University of Manchester");
  ra.delete_subkey("Test Key");

  return 0;
}


