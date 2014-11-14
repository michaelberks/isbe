#ifndef ncm_sftp_handler_h_
#define ncm_sftp_handler_h_

#include <vcl_string.h>

#include <nailfold/ncm_server_handler.h>

class ncm_sftp_handler : public ncm_server_handler
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_sftp_handler(vcl_string server = "", 
                   vcl_string username = "", 
                   vcl_string password = "");

// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
  // Members and functions visible only to objects of this class
  
  virtual vcl_string connect_string() const;

};

#endif // ncm_sftp_handler_h_