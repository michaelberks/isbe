//:
// \file
// \brief Functions for Input/Output of ncm_vessel_properties:
//          version()
//          is_a()
//          is_class()
//          vsl_b_write()
//          vsl_b_read()
//          t_write()
//          t_read()
//          operator<<
//          b_write()
//          b_read()
//          print_summary()
//          add_to_binary_loader()
// \author Phil Tresadern

#include "nailfold/ncm_vessel_properties.h"

#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_vector_io.h>
#include <vsl/vsl_indent.h>

#include <vul/vul_string.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_tuple.h>
#include <mbl/mbl_parse_sequence.h>

//
// Public functions

//
//: Version number of this class's binary IO
short ncm_vessel_properties::version() const
{ 
  return 3; 
}

//
//: The class name
vcl_string ncm_vessel_properties::is_a() const
{
  return "ncm_vessel_properties"; 
}

//
//: Class name equivalence
bool ncm_vessel_properties::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str); 
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_vessel_properties& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_vessel_properties& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_vessel_properties::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << "version: " << version() << '\n';

  // version 2
  if (is_distal_)
    tfs << vsl_indent() << "is_distal: yes" << '\n';
  else
    tfs << vsl_indent() << "is_distal: no" << '\n';

  // version 1
  tfs << vsl_indent() << "size: " << size_string() << '\n';
  tfs << vsl_indent() << "shape: " << shape_string() << '\n';

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // " << is_a() << '\n';

  tfs << vcl_flush;
}

//
//: Text read
void ncm_vessel_properties::t_read(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  clear();

  mbl_read_props_type props = mbl_read_props(tfs);

  short version = 
      vul_string_atoi(props.get_required_property("version"));

  switch (version)
  {
		case 97:
			//From auto software, fall through for now
    case 3:
    case 2:
      {
        vcl_stringstream distal_stream(props.get_required_property("is_distal"));
        is_distal_ = (distal_stream.str() == "yes");
      }
      // Fall through

    case 1:
      {
        vcl_stringstream shape_stream(props.get_required_property("shape"));
        vcl_string shape_string = shape_stream.str();
        set_shape_from_string(shape_string);

        vcl_stringstream size_stream(props.get_required_property("size"));
        vcl_string size_string = size_stream.str();

        // Relabel what were (before v3) ambiguous label combinations.
        if (version < 3)
        {
          if (is_shape_undefined() || is_shape_normal()) 
          {
            if (size_string == "Giant")
              size_string = "Undefined"; // Could be Irregular instead
          }
          else if (is_shape_tortuous() || is_shape_ramified())
          {
            if (size_string == "Giant")
              size_string = "Irregular"; // Giants can only be normal
            else if (size_string == "Enlarged")
              size_string = "Undefined"; // Could be Irregular instead
          }
        }

        set_size_from_string(size_string);

        break;
      }

    default:
      // probably a bit extreme to abort() but will leave it for now
      vcl_cerr << is_a() << "::t_read() ";
      vcl_cerr << "Unexpected version number " << version << vcl_endl;
      abort();
  }
}

//
//: Text output by reference
vcl_ostream& operator<<(vcl_ostream& os, const ncm_vessel_properties& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text output by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_vessel_properties* b)
{
  if (b)  
    return os << *b;
  else
    return os << "No " << b->is_a() << " defined.";
}

//
// Protected functions

//
//: Binary write
void ncm_vessel_properties::b_write(vsl_b_ostream& bfs) const
{
  vsl_b_write(bfs,version());

  //// version 2
  //vsl_b_write(bfs,new_variable_);

  //// version 1
  //vsl_b_write(bfs,points_);
  //vsl_b_write(bfs,widths_);
  //vsl_b_write(bfs,apex_indices_);
  //vsl_b_write(bfs,vessel_traits_);
}

//
//: Binary read
void ncm_vessel_properties::b_read(vsl_b_istream& bfs)
{
  short version;
  vsl_b_read(bfs,version);
  //switch (version)
  //{
  //  //case (2): 
  //  //  vsl_b_read(bfs,new_variable_);
  //  //  // no break

  //  case (1):
  //    vsl_b_read(bfs,points_);
  //    vsl_b_read(bfs,widths_);
  //    vsl_b_read(bfs,apex_indices_);
  //    vsl_b_read(bfs,vessel_traits_);
  //    break;

  //  default:
  //    vcl_cerr << "ncm_vessel_properties::b_read() ";
  //    vcl_cerr << "Unexpected version number " << version << vcl_endl;
  //    abort();
  //}
}

//
//: Print summary to output stream os
void ncm_vessel_properties::print_summary(vcl_ostream& os) const
{
  if (is_distal_)
    os << "Is distal: yes" << '\n';
  else
    os << "Is distal: no" << '\n';
  os << "Size: " << size_string() << '\n';
  os << "Shape: " << shape_string() << '\n';
  os << vcl_flush;
}

