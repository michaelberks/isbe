#ifndef __ncm_vessel_properties_h__
#define __ncm_vessel_properties_h__

//:
// \file
// \brief Properties of a blood vessel, including its size (normal, enlarged,
//        giant) and shape (normal, tortuous, ramified)
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_string.h>

// forward declarations
#include <vcl_iosfwd.h>
#include <vsl/vsl_fwd.h>
class ncm_vessel;

class ncm_vessel_properties
{
//  INTERFACE

public:
  //  No member variables here, please
  ncm_vessel_properties(ncm_vessel* parent = NULL);

  //: Reset to undefined properties
  void clear();

  //: Set vessel traits
  void set_distal(bool is_distal = true);

  void set_size_undefined();
  void set_size_normal();
  void set_size_enlarged();
  void set_size_giant();
  void set_size_irregular();
	void set_size_auto();

  void set_shape_undefined();
  void set_shape_normal();
  void set_shape_tortuous(bool is_tortuous = true); // i.e. nonspecific
  void set_shape_ramified(bool is_ramified = true); // i.e. angiogenic
	void set_shape_auto(); // i.e. angiogenic

  //: Get vessel traits
  bool is_distal() const;

  bool is_size_undefined() const;
  bool is_size_normal() const;
  bool is_size_enlarged() const;
  bool is_size_giant() const;
  bool is_size_irregular() const;
	bool is_size_auto() const;

  bool is_shape_undefined() const;
  bool is_shape_normal() const;
  bool is_shape_tortuous() const; // i.e. nonspecific
  bool is_shape_ramified() const; // i.e. angiogenic
	bool is_shape_auto() const; // i.e. angiogenic

  //: Copy size and shape properties from other structure
  void copy_size_from(const ncm_vessel_properties& properties);
  void copy_shape_from(const ncm_vessel_properties& properties);
  void copy_properties_from(const ncm_vessel_properties& properties);

  //: Set size and shape properties given a string
  void set_size_from_string(const vcl_string& size_string);
  void set_shape_from_string(vcl_string shape_string);

  unsigned property_value() const;
  vcl_string size_string() const;
  vcl_string shape_string() const;
  vcl_string as_string() const;

  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_vessel* clone() const;

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_vessel_properties& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_vessel_properties& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_vessel_properties& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os, const ncm_vessel_properties* b);


//  IMPLEMENTATION

protected:

  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private: // Methods
  
  unsigned vessel_size() const;
  unsigned vessel_shape() const;

  //: Tell the underlying annotation that the vessel has changed.
  void set_modified_if(bool is_modified = true);

private: // Variables

  // Vessel to which these properties apply
  ncm_vessel* parent_;

  //: Vessel descriptors
  enum vessel_descriptor {
    Vessel_Undefined = 0x0000,

    // Size
    Vessel_Normal_Size = 0x0001,  // = 1
    Vessel_Enlarged = 0x0002,     // = 2
    Vessel_Giant = 0x0004,        // = 4
    Vessel_Irregular = 0x0008,    // = 8
		Vessel_Auto = 0x0010, // = 16
    //Vessel_XXX = 0x0010, ...

    // Shape
    Vessel_Normal_Shape = 0x0100, // = 256
    Vessel_Tortuous = 0x0200,     // = 512
    Vessel_Ramified = 0x0400,     // = 1024
    //Vessel_XXX = 0x0800,
    //Vessel_XXX = 0x1000, ...

    Vessel_First = Vessel_Undefined,
    Vessel_Last = Vessel_Ramified
  }; 
  static const unsigned vessel_size_mask  = 0x00ff;
  static const unsigned vessel_shape_mask = 0xff00;

  //: Is the vessel in the distal row?
  bool is_distal_;

  //: Vessel description (OR-ed descriptors)
  unsigned vessel_properties_;
};

#endif // __ncm_vessel_properties_h__