#include "ncm_vessel_properties.h"

#include <vcl_iostream.h>

#include <mbl/mbl_parse_string_list.h>

#include <nailfold/ncm_vessel.h>

//
//  Define class member functions

//: Constructor
ncm_vessel_properties::ncm_vessel_properties(ncm_vessel* parent /* = NULL */)
: vessel_properties_(Vessel_Undefined),
  is_distal_(true),
  parent_(parent)
{
}

//
//: Reset to undefined properties
void ncm_vessel_properties::clear()
{
  set_modified_if(vessel_properties_ != Vessel_Undefined);
  vessel_properties_ = Vessel_Undefined;
}

//
//: Set vessel properties
void ncm_vessel_properties::set_distal(bool is_distal /* = true */)
{
  set_modified_if(is_distal_ != is_distal);
  is_distal_ = is_distal;
}

void ncm_vessel_properties::set_size_undefined()
{
  const unsigned new_properties = vessel_shape() | Vessel_Undefined;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_size_normal()
{
  const unsigned new_properties = vessel_shape() | Vessel_Normal_Size;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_size_enlarged()
{
  const unsigned new_properties = vessel_shape() | Vessel_Enlarged;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_size_giant()
{
  // Giant vessels have Normal shape, by definition
  const unsigned new_properties = Vessel_Normal_Shape | Vessel_Giant;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_size_irregular()
{
  const unsigned new_properties = vessel_shape() | Vessel_Irregular;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_size_auto()
{
  const unsigned new_properties = vessel_shape() | Vessel_Auto;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}

void ncm_vessel_properties::set_shape_undefined()
{
  const unsigned new_properties = vessel_size() | Vessel_Undefined;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_shape_normal()
{
  const unsigned new_properties = vessel_size() | Vessel_Normal_Shape;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}
void ncm_vessel_properties::set_shape_tortuous(bool is_tortuous /* = true */)
{
  // Giant vessels cannot be Tortuous/Nonspecific by definition.
  if (!is_size_giant())
  {
    const unsigned new_properties = vessel_size() | Vessel_Tortuous;
    set_modified_if(vessel_properties_ != new_properties);
    vessel_properties_ = new_properties;
  }
}
void ncm_vessel_properties::set_shape_ramified(bool is_ramified /* = true */)
{
  // Giant vessels cannot be Ramified/Angiogenic by definition.
  if (!is_size_giant())
  {
    const unsigned new_properties = vessel_size() | Vessel_Ramified;
    set_modified_if(vessel_properties_ != new_properties);
    vessel_properties_ = new_properties;
  }
}
void ncm_vessel_properties::set_shape_auto()
{
  const unsigned new_properties = vessel_size() | Vessel_Auto;
  set_modified_if(vessel_properties_ != new_properties);
  vessel_properties_ = new_properties;
}

//
//: Get vessel properties
bool ncm_vessel_properties::is_distal() const
{
  return is_distal_;
}

bool ncm_vessel_properties::is_size_undefined() const
{
  return (vessel_size() == Vessel_Undefined);
}
bool ncm_vessel_properties::is_size_normal() const
{
  return (vessel_size() == Vessel_Normal_Size);
}
bool ncm_vessel_properties::is_size_enlarged() const
{
  return (vessel_size() == Vessel_Enlarged);
}
bool ncm_vessel_properties::is_size_giant() const
{
  return (vessel_size() == Vessel_Giant);
}
bool ncm_vessel_properties::is_size_irregular() const
{
  return (vessel_size() == Vessel_Irregular);
}
bool ncm_vessel_properties::is_size_auto() const
{
  return (vessel_size() == Vessel_Auto);
}

bool ncm_vessel_properties::is_shape_undefined() const
{
  return ((vessel_shape() | Vessel_Undefined) == 0);
}
bool ncm_vessel_properties::is_shape_normal() const
{
  return ((vessel_shape() & Vessel_Normal_Shape) != 0);
}
bool ncm_vessel_properties::is_shape_tortuous() const
{
  return ((vessel_shape() & Vessel_Tortuous) != 0);
}
bool ncm_vessel_properties::is_shape_ramified() const
{
  return ((vessel_shape() & Vessel_Ramified) != 0);
}
bool ncm_vessel_properties::is_shape_auto() const
{
  return is_size_auto();
}
 
//: Copy size and shape properties from existing structure
void ncm_vessel_properties::copy_size_from(
  const ncm_vessel_properties& properties)
{
  if (properties.is_size_giant())
    set_size_giant(); // Also sets the shape to Normal.
  else
  {
    const unsigned new_properties = vessel_shape() | properties.vessel_size();
    set_modified_if(vessel_properties_ != new_properties);
    vessel_properties_ = new_properties;
  }
}
void ncm_vessel_properties::copy_shape_from(
  const ncm_vessel_properties& properties)
{
  if (is_size_giant())
    set_shape_normal(); // Giants are Normal by definition.
  else
  {
    const unsigned new_properties = vessel_size() | properties.vessel_shape();
    set_modified_if(vessel_properties_ != new_properties);
    vessel_properties_ = new_properties;
  }
}

//
//: Set size and shape from string input
void ncm_vessel_properties::set_size_from_string(const vcl_string& size_string)
{
  if (size_string == "Undefined")
    set_size_undefined();
  else if (size_string == "Normal")
    set_size_normal();
  else if (size_string == "Enlarged")
    set_size_enlarged();
  else if (size_string == "Giant")
    set_size_giant();
  else if (size_string == "Irregular")
    set_size_irregular();
	 else if (size_string == "auto")
    set_size_auto();
  else
  {
    // go crazy
    vcl_cerr << "Unrecognized size type: " << size_string << vcl_endl;
  }
}
void ncm_vessel_properties::set_shape_from_string(vcl_string shape_string)
{
  // Split string into its constituent parts
  // e.g. "Tortuous Ramified Bushy" -> {"Tortuous", "Ramified", "Bushy"}

  shape_string = "{ " + shape_string + " }";
  vcl_vector<vcl_string> items;
  mbl_parse_string_list(shape_string, items);

  vcl_vector<vcl_string>::const_iterator it = items.begin();
  for (it = items.begin(); it != items.end(); ++it)
  {
    if (*it == "Undefined")
    {
      set_shape_undefined();
      return; // Undefined => nothing else to say
    }
    else if (*it == "Normal")
    {
      set_shape_normal();
      return; // Normal => nothing else to say
    }
    else if ((*it == "Tortuous") || 
             (*it == "Meandering") ||
             (*it == "Non-specific"))
      set_shape_tortuous();
    else if ((*it == "Ramified") || 
             (*it == "Angiogenetic") ||
             (*it == "Angiogenic"))

      set_shape_ramified();

		else if ((*it == "auto"))
			set_shape_auto();

    else
    {
      // go crazy
      vcl_cerr << "Unrecognized shape type: " << *it << vcl_endl;
    }
  }
}

//
//: Get vessel properties as a single value
unsigned ncm_vessel_properties::property_value() const
{
  return vessel_properties_;
}

vcl_string ncm_vessel_properties::size_string() const
{
  // Vessel must be one and only one size, otherwise it is undefined

  if (is_size_normal())
    return "Normal";
  else if (is_size_enlarged())
    return "Enlarged";
  else if (is_size_giant())
    return "Giant";
  else if (is_size_irregular())
    return "Irregular";
  else
    return "Undefined";
}

vcl_string ncm_vessel_properties::shape_string() const
{
  // A vessel can have more than one shape descriptor

  if (is_shape_undefined())
    return "Undefined";

  vcl_string s = "";

  if (is_shape_normal())
    s += "Normal ";
  if (is_shape_tortuous())
    s += "Non-specific ";
  if (is_shape_ramified())
    s += "Angiogenic ";

  return s;
}

vcl_string ncm_vessel_properties::as_string() const
{
  return "Size: " + size_string() + "; Shape: " + shape_string();
}

//
//  Private methods
//

unsigned ncm_vessel_properties::vessel_size() const
{
  return (vessel_properties_ & vessel_size_mask);
}
unsigned ncm_vessel_properties::vessel_shape() const
{
  return (vessel_properties_ & vessel_shape_mask);
}
 
//: Tell the underlying annotation that the vessel has changed.
void ncm_vessel_properties::set_modified_if(bool is_modified /* = true */)
{
  if ((parent_ != NULL) && is_modified)
    parent_->set_modified();
}