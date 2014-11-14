#include "ncm_encrypt.h"

//:
// \file
// \brief A very (very) basic encryption that I use to hardcode strings with 
//        sensitive information. Once embedded in the binary code, the 
//        encrypted string shouldn't be obvious enough for anyone to recover 
//        and decrypt it.
//        It should *not* be used to encrypt anything for storing in a file!
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

vcl_string ncm_encrypt(vcl_string s)
{
  for (unsigned i = 0; i < s.length(); ++i)
  {
    if (s[i] % 2 == 0)
      s[i] += 128;
    else
      s[i] += 64;
  }
  return s;
}

vcl_string ncm_decrypt(vcl_string s)
{
  for (unsigned i = 0; i < s.length(); ++i)
  {
    if (s[i] % 2 == 0)
      s[i] -= 128;
    else
      s[i] -= 64;
  }
  return s;
}
