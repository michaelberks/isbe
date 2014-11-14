#ifndef __ncm_qimagehandler_h__
#define __ncm_qimagehandler_h__

//:
// \file
// \brief View including ability to zoom in and out
//        Adapted from Tim's qvcr_zoom_view
// \author Phil Tresadern

#include <QObject>
#include <QImage>
#include <QPixmap>
#include <QRectF>

#include <vcl_vector.h>
#include <vcl_string.h>

#include <vil/vil_image_view.h>

#include <nailfold/vil_connected_components.h>

class QGraphicsPixmapItem;

class ncm_qimagehandler : public QObject
{
  Q_OBJECT // needed if we want to use signals or slots

// INTERFACE

public:

  enum imageType { TypeUndefined = -1,
                   TypeCapillarogram = 0,
                   TypeDermatogram = 1,
  
                   TypeFirst = TypeCapillarogram,
                   TypeLast = TypeDermatogram};

  //: Default constructor
  ncm_qimagehandler(QWidget* parent = 0);

  //: Destruct
  ~ncm_qimagehandler();

  //: Clear any existing image
  void clear();

  //: Return true if the image is not currently being modified
  bool isValid() const;

  void set_filename(const vcl_string& filename);
  vcl_string filename() const;

  unsigned width() const;
  unsigned height() const;

  QPixmap const& raw_pixmap() const;
  QPixmap const& line_pixmap() const;
  QPixmap const& edge_pixmap() const;

  //: Gets nearest non-zero pixel position from supplied coordinates
  QPoint line_pixel_nearest_to(const QPointF& pos) const;
  QPoint edge_pixel_nearest_to(const QPointF& pos) const;

  bool image_is_large() const;

  void set_contrast(int min_contrast, int max_contrast,
                    bool auto_update = true);

  void set_min_contrast(int min_contrast, bool auto_update = true);
  int min_contrast() const;

  void set_max_contrast(int max_contrast, bool auto_update = true);
  int max_contrast() const;

  void set_line_size_threshold(int line_size_threshold);

  //: Get a list of all nonzero pixels in background image
  void get_line_pixels();
  void get_edge_pixels();

  imageType image_type() const;

signals:

  void imageLoaded(bool success = true);
  void rawImageChanged(bool changed = true);
  void lineImageChanged(bool changed = true);
  void edgeImageChanged(bool changed = true);

  void statusChanged(const QString& status, 
                      double progress = 0.0);

  void imageTypeChanged();

public slots:

  //: Load image specified by filename_
  void load_image();

  //: Given a rectangle, compute the min and max contrast values
  void set_contrast_from(QRect image_rect);

  //: Find ridges in the image
  void find_ridges();

  //: Find edges in the image
  void find_edges();


// IMPLEMENTATION

protected:

private:

  //: Load an image from a filename
  void load_image_from(const vcl_string& filename);
  
  //: Provide a nailfold image
  void set_image(const vil_image_view<vxl_byte>& image);

  void set_status(const vcl_string& new_status,
                  double progress = 0.0);

  //: Find 'on' pixels in binary pixmap
  void get_pixmap_pixels(const QPixmap& pixmap, 
                         vcl_vector<QPoint>& pixel_list);

  //: Find pixel nearest to given (x,y) coordinate from given list
  QPoint pixel_nearest_to(const QPointF& pos, 
                          const vcl_vector<QPoint>& pixel_list) const;

  //: Update entries of and apply colour tables
  void update_raw_colour_table();
  void apply_raw_colour_table();

  void update_centreline_colour_table();
  void apply_centreline_colour_table();
  
  void update_edge_colour_table();
  void apply_edge_colour_table();

  void set_image_type(imageType image_type);


  //  Member variables

  //: When invalid_ == 0, the image is valid and not currently undergoing any
  //  changes. For every function that writes to the class (i.e. any non-const
  //  method), invalid_ is incremented at the entry point, and decremented at
  //  every exit point. This acts like a mutex with the exception that it works
  //  over functions that are called at different levels
  int invalid_;

  //: Filename of current image
  vcl_string filename_;

  //: String indicating status of the scene
  vcl_string status_;

  //: Grey value to map to 0 (black)
  int min_contrast_;

  //: Grey value to map to 255 (white)
  int max_contrast_;

  //: Line size threshold
  int line_size_threshold_;

  //: Image type (capillarogram vs dermatogram)
  imageType image_type_;

  //: Line response threshold
  //int line_response_threshold_;

  //: Raw VXL image (after flattening and cropping)
  vil_image_view<vxl_byte> raw_vxl_image_;

  //: Raw nailfold image and its associated pixmap
  QImage raw_image_;
  QPixmap raw_pixmap_;

  //: Colormap for raw images
  QVector<QRgb> raw_colour_table_;

  //: Centreline image and its associated pixmap
  QImage line_image_;
  QPixmap line_pixmap_;

  //: Colormap for centreline components
  QVector<QRgb> centreline_colour_table_;

  //: Line image and its associated pixmap
  QImage edge_image_;
  QPixmap edge_pixmap_;

  //: Colormap for line components
  QVector<QRgb> edge_colour_table_;

  //: Class for processing connected components
  vil_connected_components conn_comp_;

  //: Nonzero pixels
  vcl_vector<QPoint> line_pixels_;
  vcl_vector<QPoint> edge_pixels_;
};

#endif // __ncm_qimagehandler_h__
