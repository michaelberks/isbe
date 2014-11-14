#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_iomanip.h>
#include <vcl_ctime.h>

#include <vil/vil_image_view.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>

#include <nailfold/ncm_image_server.h>
#include <nailfold/ncm_scp_handler.h>

//
// Main function
int main( int argc, char* argv[] )
{
  ncm_image_server server("U:/tmp/test_qmarkup/Inbox");
  vcl_cout << server.n_images() << vcl_endl;
  vcl_cout << server.n_unmarked() << vcl_endl;  

  ncm_scp_handler scp_handler;
  server.pull_images_from(scp_handler, "manchester");

  vcl_cout << server.first_image() << vcl_endl;
  bool is_valid_image = true;

  is_valid_image = true;
  while (is_valid_image)
  {
    vcl_string image_name = server.next_image();
    vcl_cout << image_name << vcl_endl;

    is_valid_image = (image_name != "");
  }

  vcl_cout << server.current_image() << vcl_endl;
  is_valid_image = true;
  while (is_valid_image)
  {
    vcl_string image_name = server.prev_image();
    vcl_cout << image_name << vcl_endl;

    is_valid_image = (image_name != "");
  }

  server.first_image();
  server.set_current_marked();
  is_valid_image = true;
  while (is_valid_image)
  {
    vcl_string image_name = server.next_unmarked_image();
    vcl_cout << image_name << vcl_endl;

    is_valid_image = (image_name != "");
  }

  server.last_image();
  is_valid_image = true;
  while (is_valid_image)
  {
    vcl_string image_name = server.prev_unmarked_image();
    vcl_cout << image_name << vcl_endl;

    is_valid_image = (image_name != "");
  }

  return 0;
}


