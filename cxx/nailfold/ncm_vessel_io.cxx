//:
// \file
// \brief Functions for Input/Output of ncm_vessel:
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

#include "ncm_vessel.h"

#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_vector_io.h>
#include <vsl/vsl_indent.h>

#include <vul/vul_string.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_tuple.h>
#include <mbl/mbl_parse_sequence.h>
#include <mbl/mbl_parse_keyword_list.h>

#include "ncm_apex.h"

//
// Public functions

//
//: Version number of this class's binary IO
short ncm_vessel::version() const
{ 
  return 1; 
}

//
//: The class name
vcl_string ncm_vessel::is_a() const
{
  return "ncm_vessel"; 
}

//
//: Class name equivalence
bool ncm_vessel::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str); 
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_vessel& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_vessel& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_vessel::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << "version: " << version() << '\n';

  //// version 2
  //tfs << new_variable_;

  // version 1
  properties().t_write(tfs);

  tfs << vsl_indent() << "anchor: { \n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << anchor().x() << " " << anchor().y() << "\n";

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // anchor\n";

  tfs << vsl_indent() << "points: {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << "// (x, y, width)\n";
  for (unsigned point = 0; point < points_.size(); ++point)
  {
    tfs << vsl_indent() << points_[point].x() << " " 
                        << points_[point].y() << " " 
                        << widths_[point] << '\n';
  }

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // points\n";

  tfs << vsl_indent() << "apices: {\n";
  vsl_indent_inc(tfs);

  for (unsigned apex = 0; apex < apices_.size(); ++apex)
    apices_[apex]->t_write(tfs);

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // apices\n";

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // " << is_a() << '\n';

  tfs << vcl_flush;
}

//
//: Text read
void ncm_vessel::t_read(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  points_.resize(0);
  widths_.resize(0);
  //apices_.resize(0)?

  mbl_read_props_type props = mbl_read_props(tfs);

  short version = 
      vul_string_atoi(props.get_required_property("version"));

  switch (version)
  {
    case 97:
			//From auto software, fall through for now..

    case 1:
      {
        vcl_stringstream props_stream(props.get_required_property("ncm_vessel_properties"));
        properties().t_read(props_stream);

        vcl_stringstream anchor_stream(props.get_required_property("anchor"));
        vcl_vector<double> anchor_data;
        mbl_parse_sequence(anchor_stream, 
                           vcl_back_inserter(anchor_data), double());
        set_anchor_to(vgl_point_2d<double>(anchor_data[0], anchor_data[1]));

        // get (x,y,width) values from points property
        vcl_stringstream points_stream(props.get_required_property("points"));

        vcl_vector<double> points_data;
        mbl_parse_sequence(points_stream, 
                           vcl_back_inserter(points_data), double());
        for (unsigned p = 2; p < points_data.size(); p+=3)
        {
          const double x = points_data[p-2];
          const double y = points_data[p-1];
          const double w = points_data[p];

          points_.push_back(vgl_point_2d<double>(x,y));
          widths_.push_back(w);
        }

        // Read apices
        vcl_string apices_data = props.get_required_property("apices");
        vcl_stringstream apices_stream(apices_data);
        vcl_vector<vcl_string> apex_strings;
        const bool discard_comments = true;
        mbl_parse_keyword_list(apices_stream, "ncm_apex:", apex_strings,
                               discard_comments);

        for (vcl_vector<vcl_string>::iterator it = apex_strings.begin();
             it != apex_strings.end(); ++it)
        {
          ncm_apex* new_apex = new ncm_apex(*this, 0, 0);
          new_apex->t_read(vcl_stringstream(*it));
          apices_.push_back(new_apex);
        }

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
vcl_ostream& operator<<(vcl_ostream& os, const ncm_vessel& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text output by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_vessel* b)
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
void ncm_vessel::b_write(vsl_b_ostream& bfs) const
{
  vsl_b_write(bfs,version());

  //// version 2
  //vsl_b_write(bfs,new_variable_);

  // version 1
  vsl_b_write(bfs,points_);
  vsl_b_write(bfs,widths_);
  //vsl_b_write(bfs,apex_indices_);
  //vsl_b_write(bfs,properties_);
}

//
//: Binary read
void ncm_vessel::b_read(vsl_b_istream& bfs)
{
  short version;
  vsl_b_read(bfs,version);
  switch (version)
  {
    //case (2): 
    //  vsl_b_read(bfs,new_variable_);
    //  // no break

    case (1):
      vsl_b_read(bfs,points_);
      vsl_b_read(bfs,widths_);
      //vsl_b_read(bfs,apex_indices_);
      //vsl_b_read(bfs,vessel_properties_);
      break;

    default:
      vcl_cerr << "ncm_vessel::b_read() ";
      vcl_cerr << "Unexpected version number " << version << vcl_endl;
      abort();
  }
}


//
//: Print summary to output stream os
void ncm_vessel::print_summary(vcl_ostream& os) const
{
  os << "Points: " << points_.size() << '\n';
  //os << "Apices: " << apex_indices_.size() << '\n';
  os << "Traits: ";
  if (properties().is_size_enlarged()) 
    os << "Enlarged ";
  if (properties().is_size_giant()) 
    os << "Giant ";
  if (properties().is_shape_tortuous()) 
    os << "Tortuous ";
  os << '\n';

  os << vcl_flush;
}

