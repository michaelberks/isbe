#ifndef ncm_registry_accessor_h_
#define ncm_registry_accessor_h_

#include <vcl_string.h>

#include <windows.h>

class ncm_registry_accessor
{
public:
  // Public member functions

  //: Default constructor
  ncm_registry_accessor();

  //: Destructor
  ~ncm_registry_accessor();

  //: Open HKEY_CURRENT_USER key.
  long get_current_user_key();

  //: Open a subkey of the current subkey.
  long get_subkey(vcl_string subkey);

  //: Create a new subkey of the current subkey.
  long create_subkey(vcl_string subkey);

  //: Delete a new subkey of the current subkey.
  long delete_subkey(vcl_string subkey);

  //: Write or create a string name/value pair within the current subkey.
  long write_string(vcl_string name, vcl_string value);

  //: Read a string name/value pair within the current subkey.
  long read_string(vcl_string name, vcl_string& value);

  //: Write or create a numeric name/value pair within the current subkey.
  long write_numeric(vcl_string name, int value);

  //: Read a numeric name/value pair within the current subkey.
  long read_numeric(vcl_string name, int& value);

  //: Write or create a numeric name/value pair within the current subkey.
  long write_boolean(vcl_string name, bool value);

  //: Read a numeric name/value pair within the current subkey.
  long read_boolean(vcl_string name, bool& value);

protected:
  // Protected member functions

private:
  // Private member functions

  // Private member data

  //: Handle to the currently open subkey.
  HKEY handle_;
};

#endif // ncm_registry_accessor_h_