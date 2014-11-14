#include <vcl_iostream.h>

#include <vcl_vector.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>

// forward declaration
class myclass
{
public:
  myclass() : ref(10), p_ref(&ref) {}
  int* myfunc() { return p_ref; }
  int* const c_myfunc() { return p_ref; }

private:
  int ref;
  int* p_ref;
};


int ncm_vessel_test()
{
  myclass c;
  int* const p = c.myfunc();

  //int*& p1 = c.myfunc();
  //int*& p2 = c.c_myfunc();
  int* const& p3 = c.myfunc();
  int* const& p4 = c.c_myfunc();

  vcl_cout << *p3 << vcl_endl;

  return 0;
  ncm_annotation annotation;
  vgl_point_2d<double> anchor(1.0, 2.0);

  ncm_vessel* v = annotation.create_vessel_at(anchor.x(), anchor.y());

  if (v->anchor() != anchor)
    return 1;

  for (unsigned i = 0; i < 5; ++i)
    v->add_venous_point(0, i);

  for (unsigned i = 0; i < v->n_points(); ++i)
    vcl_cout << *(v->point(i)) << vcl_endl;

  vgl_vector_2d<double> n = v->normal_at(0,2);
  vcl_cout << n << vcl_endl;

  annotation.vessel(0)->add_arterial_point(0,0);

  const ncm_annotation& alias = annotation;
//  alias.vessel(0)->add_arterial_point(0,0);

  const vgl_point_2d<double> * epp = v->venous_endpoint();
  epp = v->venous_endpoint();

  for (unsigned i = 0; i < v->n_points(); ++i)
    vcl_cout << *(v->point(i)) << vcl_endl;

  return 0;
}

int main()
{
  //ncm_annotation_test();
  ncm_vessel_test();

  return 0;
}

