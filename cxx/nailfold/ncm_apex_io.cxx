//:
// \file
// \brief Functions for Input/Output of ncm_apex:
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

#include "nailfold/ncm_apex.h"

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
short ncm_apex::version() const
{ 
  return 1; 
}

//
//: The class name
vcl_string ncm_apex::is_a() const
{
  return "ncm_apex"; 
}

//
//: Class name equivalence
bool ncm_apex::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str); 
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_apex& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_apex& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_apex::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << "version: " << version() << '\n';

  //// version 2
  //tfs << new_variable_;

  // version 1
  tfs << vsl_indent() << "inner_point: {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << inner_point().x() << " " 
                      << inner_point().y() << "\n";
  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // inner_point\n";

  tfs << vsl_indent() << "outer_point: {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << outer_point().x() << " " 
                      << outer_point().y() << " \n";
  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // outer_point\n";

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // " << is_a() << '\n';

  tfs << vcl_flush;
}

//
//: Text read
void ncm_apex::t_read(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  mbl_read_props_type props = mbl_read_props(tfs);

  short version = 
      vul_string_atoi(props.get_required_property("version"));

  switch (version)
  {
    case 97:
			//From auto software, fall through for now
    case 1:
      {
        vcl_stringstream inner_stream(props.get_required_property("inner_point"));
        vcl_vector<double> inner_data;
        mbl_parse_sequence(inner_stream, 
                           vcl_back_inserter(inner_data), double());
        set_inner_point(inner_data[0], inner_data[1]);

        vcl_stringstream outer_stream(props.get_required_property("outer_point"));
        vcl_vector<double> outer_data;
        mbl_parse_sequence(outer_stream, 
                           vcl_back_inserter(outer_data), double());
        set_outer_point(outer_data[0], outer_data[1]);
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
vcl_ostream& operator<<(vcl_ostream& os, const ncm_apex& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text output by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_apex* b)
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
void ncm_apex::b_write(vsl_b_ostream& bfs) const
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
void ncm_apex::b_read(vsl_b_istream& bfs)
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
  //    vcl_cerr << "ncm_apex::b_read() ";
  //    vcl_cerr << "Unexpected version number " << version << vcl_endl;
  //    abort();
  //}
}

//
//: Print summary to output stream os
void ncm_apex::print_summary(vcl_ostream& os) const
{
  //os << "Points: " << points_.size() << '\n';
  //os << "Apices: " << apex_indices_.size() << '\n';
  //os << "Traits: ";
  //if (is_enlarged()) 
  //  os << "Enlarged ";
  //if (is_giant()) 
  //  os << "Giant ";
  //if (is_bushy()) 
  //  os << "Bushy ";
  //os << '\n';

  //os << vcl_flush;
}

