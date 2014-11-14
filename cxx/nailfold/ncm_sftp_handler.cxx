#include "ncm_sftp_handler.h"

#include <nailfold/ncm_encrypt.h>

//
// Public member functions
//

ncm_sftp_handler::ncm_sftp_handler(vcl_string server /* = "" */, 
                                   vcl_string username /* = "" */, 
                                   vcl_string password /* = "" */)
: ncm_server_handler(server, username, password)
{
  connect();
}

//
// Private member functions
//

vcl_string ncm_sftp_handler::connect_string() const
{
  return ncm_decrypt("ð³£ð m³æôð mÐ ²² mð· " + password_);
}

