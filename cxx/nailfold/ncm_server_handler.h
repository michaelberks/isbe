#ifndef ncm_server_handler_h_
#define ncm_server_handler_h_

#include <vcl_string.h>

class ncm_server_handler
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_server_handler(vcl_string server = "", 
                     vcl_string username = "", 
                     vcl_string password = "");

  void set_server(vcl_string server);
  void set_username(vcl_string username);
  void set_password(vcl_string password);

  bool connect() const;
  bool is_connected() const;

  bool push(vcl_string src, vcl_string dest, vcl_string options = "") const;
  bool pull(vcl_string src, vcl_string dest, vcl_string options = "") const;


// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

  vcl_string password_;

private:
  // Members and functions visible only to objects of this class
  
  virtual vcl_string connect_string() const = 0;
  vcl_string server_string() const;

  vcl_string server_;
  vcl_string username_;
};

#endif // ncm_server_handler_h_