#include "ncm_server_handler.h"

#include <nailfold/ncm_encrypt.h>

//
// Public member functions
//

ncm_server_handler::ncm_server_handler(vcl_string server /* = "" */, 
                                       vcl_string username /* = "" */, 
                                       vcl_string password /* = "" */)
: server_(server),
  username_(username),
  password_(password)
{
}

void ncm_server_handler::set_server(vcl_string server)
{
  server_ = server;
}

void ncm_server_handler::set_username(vcl_string username)
{
  username_ = username;
}

void ncm_server_handler::set_password(vcl_string password)
{
  password_ = password;
}

// 
//: Connect to the server for the first time, automatically storing the host
//  key in the registry.
bool ncm_server_handler::connect() const
{
  // Return if any of the details are not given
  if (server_ == "" || username_ == "" || password_ == "")
    return false;

  // First, create a text file with the character 'y' and nothing else.
  // This will store the host key in the registry.
  bool success = (system("echo y > accept.txt") == 0);

  if (success)
  {
    // accept.txt created successfully.
    // Now pull connection_test.txt from server.
    vcl_string commandString = connect_string() + " " +
                               server_string() + 
                               "\"connection_test.txt\" " +
                               "\"tmp.txt\" < accept.txt";
    success = (system(commandString.c_str()) == 0);
  }

  if (success)
  {
    // connection_test.txt read successfully.
    // Now clean up.
    system("del accept.txt");
    system("del tmp.txt");
  }

  return success;
}

bool ncm_server_handler::is_connected() const
{
  bool success = pull("connection_test.txt", "./tmp.txt");

  if (success)
  {
    system("del tmp.txt");
    return true;
  }
  else
    return false;
}

bool ncm_server_handler::push(vcl_string src, vcl_string dest,
                              vcl_string options /* = "" */) const
{
  vcl_string commandString = connect_string() + " " +
                             ncm_decrypt(options) + " " + 
                             "\"" + src + "\" " +
                             server_string() + "\"" + dest + "\" ";

  return (system(commandString.c_str()) == 0);
}

bool ncm_server_handler::pull(vcl_string src, vcl_string dest, 
                              vcl_string options /* = "" */) const
{
  vcl_string commandString = connect_string() + " " +
                             ncm_decrypt(options) + " " + 
                             server_string() + "\"" + src + "\" " +
                             "\"" + dest + "\" ";

  return (system(commandString.c_str()) == 0);
}

//
// Private member functions
//

//: Return string in form <username>@<server>:
vcl_string ncm_server_handler::server_string() const
{
  return ncm_decrypt(username_ + "À" + server_ + "º");
}

