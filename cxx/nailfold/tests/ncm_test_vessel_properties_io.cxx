#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_string.h>

#include <nailfold/ncm_vessel_properties.h>

//
// Main function
int main( int argc, char* argv[] )
{
  ncm_vessel_properties vp;

  vcl_vector<vcl_string> sizes;
  sizes.push_back("Undefined");
  sizes.push_back("Normal");
  sizes.push_back("Enlarged");
  sizes.push_back("Giant");
  sizes.push_back("Irregular");

  vcl_vector<vcl_string> shapes;
  shapes.push_back("Undefined");
  shapes.push_back("Normal");
  shapes.push_back("Non-specific");
  shapes.push_back("Angiogenic");

  vcl_vector<vcl_string> yes_no;
  yes_no.push_back("yes");
  yes_no.push_back("no");
  
  for (unsigned version = 2; version <= 3; ++version)
  {
    vcl_cout << "VERSION = " << version << vcl_endl;

    for (unsigned idistal = 0; idistal < 2; ++idistal)
    {
      vcl_cout << "DISTAL = " << yes_no[idistal] << vcl_endl;

      for (unsigned isize = 0; isize < sizes.size(); ++isize)
      {
        for (unsigned ishape = 0; ishape < shapes.size(); ++ishape)
        {
          vcl_stringstream ss;
          ss<< "{" << vcl_endl
            << "  version: " << version << vcl_endl
            << "  is_distal: " << yes_no[idistal] << vcl_endl
            << "  size: " << sizes[isize] << vcl_endl
            << "  shape: " << shapes[ishape] << vcl_endl
            << "}" << vcl_endl;

          ncm_vessel_properties ivp;
          ivp.t_read(ss);

          vcl_cout << sizes[isize] << "+" << shapes[ishape] 
                   << " -> "
                   << ivp.size_string() << "+" << ivp.shape_string() 
                   << vcl_endl;
        }
      }
    }
  }

  return 0;
}