//:
// \file
// \brief Functions for Input/Output of ncm_annotation:
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

#include "nailfold/ncm_annotation.h"

#include <vcl_sstream.h>
#include <vcl_iomanip.h>

#include <vsl/vsl_binary_io.h>
//#include <vsl/vsl_binary_loader.h>
#include <vsl/vsl_indent.h>

#include <vul/vul_string.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_keyword_list.h>

#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_haemorrhage.h>

//
// Public functions

//
//: Version number of this class's binary IO
short ncm_annotation::version() const
{ 
  return version_; //
}
short ncm_annotation::latest_manual_version() const
{ 
  return 2; // QMarkup v1.4
}
short ncm_annotation::latest_auto_version() const
{ 
  return 97; // Auto markup from Matlab software
}

//
//: The class name
vcl_string ncm_annotation::is_a() const
{
  return "ncm_annotation"; 
}

//
//: Class name equivalence
bool ncm_annotation::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str);
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_annotation& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_annotation& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_annotation::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {" << '\n';

  vsl_indent_inc(tfs);
  tfs << vsl_indent() << "version: " << version() << '\n';

  //// version 3
  //vsl_b_write(bfs,new_variable_);

  // version 2
    // No change - just used as a marker to keep track of the qmarkup 
    // software being used.

  // version 1
  tfs << vsl_indent() << "observer: " << observer_name_ << '\n';
  tfs << vsl_indent() << "timestamp: " << date_string(); // EOL appended already
  tfs << vsl_indent() << "time_taken: " << time_taken_ << '\n';
  tfs << vsl_indent() << "image_grade: " << image_grade().as_string() << '\n';

  tfs << vsl_indent() << "vessels: {" << '\n';
  vsl_indent_inc(tfs);

  for (vcl_vector<ncm_vessel*>::const_iterator it = vessels_.begin();
       it != vessels_.end(); ++it)
    (*it)->t_write(tfs);
  
  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // vessels" << '\n';

  tfs << vsl_indent() << "haemorrhages: {" << '\n';
  vsl_indent_inc(tfs);

  for (vcl_vector<ncm_haemorrhage*>::const_iterator it = haemorrhages_.begin();
       it != haemorrhages_.end(); ++it)
    (*it)->t_write(tfs);
  
  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // haemorrhages" << '\n';

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // " << is_a() << '\n';
  tfs << vcl_flush;
}

//
//: Text read (from stream, without identifying token)
void ncm_annotation::t_read(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  // Disable snapshots while we load in new data
  set_snapshots(false);
  clear();

  mbl_read_props_type props = mbl_read_props(tfs);

  version_ = 
      vul_string_atoi(props.get_required_property("version"));

  switch (version_)
  {
		case 97:
			//This is the automated software (97 being the ascii conversion of "a")
			//For now try falling through and seeing what the manual reader does...
    case 2:
      // No change - just used as a marker to keep track of the qmarkup 
      // software being used.

    case 1:
      {
        observer_name_ =  props.get_required_property("observer");

        // Timestamp is not read here as it will be reset to reflect the last
        // time that an image was edited

        time_taken_ = vul_string_atof(props.get_required_property("time_taken"));

        const vcl_string grade_string = 
            props.get_required_property("image_grade");
        image_grade().set_from_string(grade_string);

        // Read vessels
        vcl_string vessels_data = props.get_required_property("vessels");
        vcl_stringstream vessels_stream(vessels_data);
        vcl_vector<vcl_string> vessel_strings;
        const bool discard_comments = true;
        mbl_parse_keyword_list(vessels_stream, "ncm_vessel:", vessel_strings,
                               discard_comments);

        for (vcl_vector<vcl_string>::iterator it = vessel_strings.begin();
             it != vessel_strings.end(); ++it)
        {
          ncm_vessel* new_vessel = ncm_vessel_handler::new_vessel(*this, 0, 0);
          new_vessel->t_read(vcl_stringstream(*it));
          vessels_.push_back(new_vessel);
        }

        // Read haemorrhages
        vcl_string haemorrhages_data = props.get_required_property("haemorrhages");
        vcl_stringstream haemorrhages_stream(haemorrhages_data);
        vcl_vector<vcl_string> haemorrhage_strings;
        mbl_parse_keyword_list(haemorrhages_stream, "ncm_haemorrhage:", haemorrhage_strings,
                               discard_comments);

        for (vcl_vector<vcl_string>::iterator it = haemorrhage_strings.begin();
             it != haemorrhage_strings.end(); ++it)
        {
          ncm_haemorrhage* new_haemorrhage = new ncm_haemorrhage(*this, 0, 0);
          new_haemorrhage->t_read(vcl_stringstream(*it));
          haemorrhages_.push_back(new_haemorrhage);
        }

        break;
      }

    default:
      // probably a bit extreme to abort() but will leave it for now
      vcl_cerr << is_a() << "::t_read() ";
      vcl_cerr << "Unexpected version number " << version_ << vcl_endl;
      abort();
  }

  update_indices();

  is_modified_ = false;

  set_snapshots(true);
}

//
//: Text read (from file, with identifying token)
void ncm_annotation::t_read(const vcl_string& filename)
{
  vcl_ifstream tfs(filename.c_str());

  mbl_read_props_type props = mbl_read_props(tfs);
  vcl_string token = is_a();
  vcl_string markup_string = props.get_required_property(token);
  vcl_stringstream markup_stream(markup_string);
  t_read(markup_stream);

  tfs.close();
}

//
//: Text write by reference
vcl_ostream& operator<<(vcl_ostream& os, const ncm_annotation& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text write by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_annotation* b)
{
  if (b)
    return os << *b;
  else
    return os << "No " << b->is_a() << " defined.";
}


//
// Protected functions

//
//: Print summary to output stream os
void ncm_annotation::print_summary(vcl_ostream& os) const
{
  os << "Observer: " << observer_name_ << '\n';
  os << "Date: " << "<date>" << '\n';
  os << "Time taken: " << time_taken_ << '\n';
  os << "Image grade: " << image_grade().as_string() << '\n';
  os << "nVessels: " << vessels_.size() << '\n';
  os << vcl_flush;
}

//
//: Binary write
void ncm_annotation::b_write(vsl_b_ostream& bfs) const
{
  vsl_b_write(bfs,version());

  //// version 2
  //vsl_b_write(bfs,new_variable_);

  // version 1
  vsl_b_write(bfs,observer_name_);
  vsl_b_write(bfs,datetime_);
  vsl_b_write(bfs,time_taken_);
  vsl_b_write(bfs,vessels_.size());
  for (vcl_vector<ncm_vessel*>::const_iterator it = vessels_.begin();
       it != vessels_.end(); ++it)
    vsl_b_write(bfs,**it);
  //vsl_b_write(bfs,haemorrhages_);
}

//
//: Binary read
void ncm_annotation::b_read(vsl_b_istream& bfs)
{
  short version;
  vsl_b_read(bfs,version);
	version_ = version;
  switch (version)
  {
		case (97):

    case (2): 
    //  vsl_b_read(bfs,new_variable_);
    //  // no break

    case (1):
      vsl_b_read(bfs,observer_name_);
      vsl_b_read(bfs,datetime_);
      vsl_b_read(bfs,time_taken_);

      unsigned n_vessels;
      vsl_b_read(bfs,n_vessels);
      vessels_.resize(n_vessels);
      for (unsigned i = 0; i < n_vessels; ++i)
      {
        vessels_[i] = ncm_vessel_handler::new_vessel(*this, 0, 0);
        vsl_b_read(bfs,*vessels_[i]);
      }

      unsigned n_haemorrhages;
      vsl_b_read(bfs,n_haemorrhages);
      vessels_.resize(n_haemorrhages);
      for (unsigned i = 0; i < n_haemorrhages; ++i)
      {
        //haemorrhages_[i] = new ncm_haemorrhage(*this, 0, 0);
        //vsl_b_read(bfs,*vessels_[i]);
      }
      break;

    default:
      vcl_cerr << is_a() << "::b_read() ";
      vcl_cerr << "Unexpected version number " << version << vcl_endl;
      abort();
  }
}

//
//: Add to list of binary loaders (for polymorphic IO)
//void vsl_add_to_binary_loader(const ncm_annotation& b)
//{
//  vsl_binary_loader<ncm_annotation>::instance().add(b);
//}
