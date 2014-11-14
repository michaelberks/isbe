#include "ncm_registry_accessor.h"

#include <vcl_iostream.h>

//:
// \file
// \brief A class that provides easy access to Windows registry entries.
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

//
// Public member functions
//

ncm_registry_accessor::ncm_registry_accessor()
{
}

ncm_registry_accessor::~ncm_registry_accessor()
{
  if (RegCloseKey(handle_) != ERROR_SUCCESS)
  {
    // handle the error?
  }
}

//: Read from HKEY_CURRENT_USER
long ncm_registry_accessor::get_current_user_key()
{
  return RegOpenCurrentUser(KEY_READ, &handle_);
}

long ncm_registry_accessor::get_subkey(vcl_string subkey)
{
  HKEY new_handle = NULL;
  LONG result = RegOpenKeyEx(handle_,
                             subkey.c_str(),
                             NULL, // Reserved. Must be zero.
                             KEY_READ,
                             &new_handle);

  if (result != ERROR_SUCCESS)
    return result;

  handle_ = new_handle;

  return ERROR_SUCCESS;
}

long ncm_registry_accessor::create_subkey(vcl_string subkey)
{
  LONG result = RegCreateKeyEx(handle_,
                               subkey.c_str(),
                               0, // Reserved. Must be zero.
                               NULL, // Class type (not used).
                               REG_OPTION_NON_VOLATILE,
                               KEY_WRITE,
                               NULL, // not used
                               &handle_,
                               NULL); // not used

  if (result != ERROR_SUCCESS)
    return result;

  return ERROR_SUCCESS;
}

long ncm_registry_accessor::delete_subkey(vcl_string subkey)
{
  // RegDeleteKeyEx() is only available under 64-bit XP
  //LONG result = RegDeleteKeyEx(handle_,
  //                             subkey.c_str(),
  //                             KEY_WOW64_32KEY,
  //                             0); // Reserved. Must be zero.

  LONG result = RegDeleteKey(handle_,
                             subkey.c_str());

  if (result != ERROR_SUCCESS)
    return result;

  return ERROR_SUCCESS;
}

//
//: Write or create a string name/value pair within the current subkey.
long ncm_registry_accessor::write_string(vcl_string name, vcl_string value)
{
  LONG result = ERROR_SUCCESS;

  // Open key with write permissions
  result = RegOpenKeyEx(handle_,
                        "", // Use existing handle
                        NULL, // Reserved. Must be zero.
                        KEY_WRITE,
                        &handle_);

  if (result != ERROR_SUCCESS)
    return result;

  // Write the value
  result = RegSetValueEx(handle_,
                         name.c_str(),
                         0, // Reserved. Must be zero.
                         REG_SZ,
                         (LPBYTE) value.c_str(), // Null-terminated version
                         value.size()+1); // Size, including terminator

  if (result != ERROR_SUCCESS)
    return result;

  return ERROR_SUCCESS;
}

//
//: Read a string name/value pair within the current subkey.
long ncm_registry_accessor::read_string(vcl_string name, vcl_string& value)
{
  LONG result = ERROR_SUCCESS;

  // Open key with read permissions
  result = RegOpenKeyEx(handle_,
                        "", // Use existing handle
                        NULL, // Reserved. Must be zero.
                        KEY_READ,
                        &handle_);

  if (result != ERROR_SUCCESS)
    return result;

  // Create buffer for retrieved data.
  DWORD dataSize = 256;
  TCHAR* pData = new TCHAR [dataSize];
  DWORD dataType;

  // Retrieve the data.
  result = RegQueryValueEx(handle_,
                           name.c_str(),
                           NULL, // Reserved. Must be NULL.
                           &dataType, // Data type - not required here.
                           (LPBYTE) pData,
                           &dataSize);

  if (dataType != REG_SZ)
  {
    // Throw an error?
  }

  if (result != ERROR_SUCCESS)
  {
    delete [] pData;
    return result;
  }

  // Copy retrieved data to value (minus the terminator) and clean up.
  value = vcl_string(pData, dataSize-1);
  delete [] pData;

  return ERROR_SUCCESS;
}

//
//: Write or create a numeric name/value pair within the current subkey.
long ncm_registry_accessor::write_numeric(vcl_string name, int value)
{
  LONG result = ERROR_SUCCESS;

  // Open key with write permissions
  result = RegOpenKeyEx(handle_,
                        "", // Use existing handle
                        NULL, // Reserved. Must be zero.
                        KEY_WRITE,
                        &handle_);

  if (result != ERROR_SUCCESS)
    return result;

  // Write the value
  result = RegSetValueEx(handle_,
                         name.c_str(),
                         0, // Reserved. Must be zero.
                         REG_DWORD,
                         (LPBYTE) &value,
                         4); // DWORD is 4 bytes

  if (result != ERROR_SUCCESS)
    return result;

  return ERROR_SUCCESS;
}

//
//: Read a numeric name/value pair within the current subkey.
long ncm_registry_accessor::read_numeric(vcl_string name, int& value)
{
  LONG result = ERROR_SUCCESS;

  // Open key with read permissions
  result = RegOpenKeyEx(handle_,
                        "", // Use existing handle
                        NULL, // Reserved. Must be zero.
                        KEY_READ,
                        &handle_);

  if (result != ERROR_SUCCESS)
    return result;

  // Create buffer for retrieved data.
  DWORD dataSize = 256;
  DWORD dataType;

  // Retrieve the data.
  result = RegQueryValueEx(handle_,
                           name.c_str(),
                           NULL, // Reserved. Must be NULL.
                           &dataType, // Data type - not required here.
                           (LPBYTE) &value,
                           &dataSize);

  if (result != ERROR_SUCCESS)
    return result;

  if (dataType != REG_DWORD)
  {
    // Throw an error?
  }

  return ERROR_SUCCESS;
}


//
//: Write or create a bool name/value pair within the current subkey.
long ncm_registry_accessor::write_boolean(vcl_string name, bool value)
{
  LONG result = ERROR_SUCCESS;

  // Open key with write permissions
  result = RegOpenKeyEx(handle_,
                        "", // Use existing handle
                        NULL, // Reserved. Must be zero.
                        KEY_WRITE,
                        &handle_);

  if (result != ERROR_SUCCESS)
    return result;

  DWORD dword_value = 0;

  if (value)
    dword_value = 1;

  // Write the value
  result = RegSetValueEx(handle_,
                         name.c_str(),
                         0, // Reserved. Must be zero.
                         REG_DWORD,
                         (LPBYTE) &dword_value,
                         4); // DWORD is 4 bytes

  if (result != ERROR_SUCCESS)
    return result;

  return ERROR_SUCCESS;
}

//
//: Read a bool name/value pair within the current subkey.
long ncm_registry_accessor::read_boolean(vcl_string name, bool& value)
{
  LONG result = ERROR_SUCCESS;

  // Open key with read permissions
  result = RegOpenKeyEx(handle_,
                        "", // Use existing handle
                        NULL, // Reserved. Must be zero.
                        KEY_READ,
                        &handle_);

  if (result != ERROR_SUCCESS)
    return result;

  // Create buffer for retrieved data.
  DWORD dataSize = 256;
  DWORD dataType;
  DWORD dword_value;

  // Retrieve the data.
  result = RegQueryValueEx(handle_,
                           name.c_str(),
                           NULL, // Reserved. Must be NULL.
                           &dataType, // Data type - not required here.
                           (LPBYTE) &dword_value,
                           &dataSize);

  if (result != ERROR_SUCCESS)
    return result;

  if (dataType != REG_DWORD)
  {
    // Throw an error?
  }

  value = (dword_value != 0);

  return ERROR_SUCCESS;
}
