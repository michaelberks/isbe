#include "ncm_qseries_options.h"

// VXL libraries
#include <vcl_iostream.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>
#include <nailfold/qtools/ncm_qvesselitem.h>

// Qt libraries
#include <QColorDialog>

//
//: Constructor
ncm_qseries_options::ncm_qseries_options(ncm_qseries_preferences& preferences,
                                         QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags),
  current_preferences_(preferences)
{
  // setup the UI
  ui.setupUi(this);

  // Take a copy of the incoming preferences that we can modify
  new_preferences_ = preferences;

  // Create a button group for the user level (Basic/Advanced)
  userLevelGroup_ = new QButtonGroup(this);
  userLevelGroup_->addButton(ui.basicUserRadio, 0);
  userLevelGroup_->addButton(ui.advancedUserRadio, 1);

  ui.gridPixelsEdit->setValidator(&doubleValidator_);
  ui.gridMillimetresEdit->setValidator(&doubleValidator_);

  // update values from environment
  read_values();
}

//
//: Destructor
ncm_qseries_options::~ncm_qseries_options()
{
}

//
//  Public slots
//

//
//: What to do when the user hits 'OK'
void ncm_qseries_options::accept()
{
  // save new settings
  write_values();

  // close dialog with return code 0 (success)
  done(0);
}

//
//  Private slots
//

//
//: Colour selection push buttons
void ncm_qseries_options::on_apexColorButton_clicked()
{
  const QRgb new_rgb = QColorDialog::getColor().rgb();
  apply_colour_to(ui.apexColorButton, new_rgb);
}
void ncm_qseries_options::on_vessel_colour_clicked()
{
  const QRgb new_rgb = QColorDialog::getColor().rgb();
  apply_colour_to(ui.vessel_colour, new_rgb);
}
void ncm_qseries_options::on_vessel_colour_selected_clicked()
{
  const QRgb new_rgb = QColorDialog::getColor().rgb();
  apply_colour_to(ui.vessel_colour_selected, new_rgb);
}

//: Grid spacing definitions
void ncm_qseries_options::on_gridPixelsEdit_textEdited(const QString& text)
{
  updateGridLabels();
}
void ncm_qseries_options::on_gridMillimetresEdit_textEdited(const QString& text)
{
  updateGridLabels();
}

//
//  Private functions
//

//
//: Give pushbutton a block of specified colour
void ncm_qseries_options::apply_colour_to(QPushButton* const button, 
                                          QRgb const& rgb)
{
  QPixmap pixmap(16,16);
  pixmap.fill(QColor(rgb));
  button->setIcon(QIcon(pixmap));
}

//
//: Get colour currently assigned to a pushbutton
QRgb ncm_qseries_options::get_colour_from(QPushButton const* const button)
{
  QPixmap pixmap = button->icon().pixmap(16,16);
  return pixmap.toImage().pixel(0,0);
}

void ncm_qseries_options::updateGridLabels()
{
  const double gridPixels = ui.gridPixelsEdit->text().toDouble();
  const double gridMillimetres = ui.gridMillimetresEdit->text().toDouble();

  const double gridPixelsPerMm = gridPixels / gridMillimetres;
  const double gridMmPerPixel = gridMillimetres / gridPixels;

  ui.gridPixelsLabel->setText(QString::number(gridPixelsPerMm));
  ui.gridMillimetresLabel->setText(QString::number(gridMmPerPixel));
}

//
//: Update dialog controls based on given preferences
void ncm_qseries_options::read_values()
{
  // alias for convenience
  ncm_qseries_preferences& cp = current_preferences_;
  
  // User properties
  ui.usernameEdit->setText(cp.username().c_str());
  if (cp.is_advanced_user())
    ui.advancedUserRadio->setChecked(true);
  else
    ui.basicUserRadio->setChecked(true);

  // Vessel appearance properties
  ncm_qvesselitem_appearance const& v_app = QGraphicsVesselItem::appearance();
  ui.placeholder_radius->setValue(v_app.placeholder_radius() * 1000);
  ui.line_width->setValue(v_app.line_width() * 1000);
  ui.scale_enlarged->setValue(v_app.scale_enlarged());
  ui.scale_giant->setValue(v_app.scale_giant());
  apply_colour_to(ui.vessel_colour, v_app.colour_normal());
  apply_colour_to(ui.vessel_colour_selected, v_app.colour_selected());
  ui.vessel_opacity->setValue(v_app.opacity() * 100.0 / 255);
  ui.vessel_opacity_selected->setValue(v_app.opacity_selected() * 100.0 / 255);
  ui.vessel_outline_only->setChecked(v_app.outline_when_selected());

  // Vessel path smoothing properties
  ui.minimum_interpoint_distance->setValue(ncm_vessel::minimum_inter_point_distance());
  ui.maximum_interpoint_distance->setValue(ncm_vessel::maximum_inter_point_distance());

  // Sorting by x-coordinate
  ui.sortVessels->setChecked(ncm_annotation::vessels_are_sorted());

  // Zoom options
  ui.labelImageZoom->setCurrentIndex(static_cast<int>(cp.label_image_zoom_));
  ui.addVesselsZoom->setCurrentIndex(static_cast<int>(cp.add_vessels_zoom_));
  ui.setVesselPropsZoom->setCurrentIndex(static_cast<int>(cp.set_vessel_props_zoom_));
  ui.addHaemorrhagesZoom->setCurrentIndex(static_cast<int>(cp.add_haemorrhages_zoom_));

  // Apex labelling zoom factors
  ui.normalWindowSize->setValue(cp.normal_vessel_zoom());
  ui.enlargedZoomRel->setValue(cp.enlarged_zoom_relative());
  ui.giantZoomRel->setValue(cp.giant_zoom_relative());

  // Vessel drawing y-offset
  ui.pathVerticalShift->setValue(cp.vessel_path_yshift() * 100);

  // Get grid spacing properties
  ui.gridSpacingSpin->setValue(cp.grid_spacing_mm());
  ui.gridPixelsEdit->setText(QString::number(cp.grid_pixels_per_mm()));
  ui.gridMillimetresEdit->setText("1.0");
  updateGridLabels();

  // Warnings
  ui.warnUngraded->setChecked(cp.warn_if_ungraded);
  ui.warnVesselsMissing->setChecked(cp.warn_if_vessels_missing);
  ui.warnSizesMissing->setChecked(cp.warn_if_sizes_missing);
  ui.warnSizeOverwrite->setChecked(v_app.warn_on_overwrite_size);
  ui.warnShapesMissing->setChecked(cp.warn_if_shapes_missing);
  ui.warnShapeOverwrite->setChecked(v_app.warn_on_overwrite_shape);
  ui.warnApicesMissing->setChecked(cp.warn_if_apices_missing);
  ui.warnPathsMissing->setChecked(cp.warn_if_paths_missing);
  ui.warnHaemorrhagesMissing->setChecked(cp.warn_if_haemorrhages_missing);
}

//
//: Get existing values from environment
void ncm_qseries_options::write_values()
{
  // aliases for convenience
  ncm_qseries_preferences& cp = current_preferences_;

  // User properties
  cp.set_username(ui.usernameEdit->text().toStdString());
  cp.set_advanced_user(ui.advancedUserRadio->isChecked());

  // Vessel appearance properties
  ncm_qvesselitem_appearance& v_app = QGraphicsVesselItem::appearance();
  v_app.set_placeholder_radius(ui.placeholder_radius->value() * 0.001);
  v_app.set_line_width(ui.line_width->value() * 0.001);
  v_app.set_scale_enlarged(ui.scale_enlarged->value());
  v_app.set_scale_giant(ui.scale_giant->value());
  v_app.set_colour_normal(get_colour_from(ui.vessel_colour));
  v_app.set_colour_selected(get_colour_from(ui.vessel_colour_selected));
  v_app.set_opacity(ui.vessel_opacity->value() * 0.01 * 255);
  v_app.set_opacity_selected(ui.vessel_opacity_selected->value() * 0.01 * 255);
  v_app.set_outline_when_selected(ui.vessel_outline_only->isChecked());

  // Vessel path smoothing properties
  ncm_vessel::set_minimum_inter_point_distance(ui.minimum_interpoint_distance->value());
  ncm_vessel::set_maximum_inter_point_distance(ui.maximum_interpoint_distance->value());

  // Sorting by x-coordinate
  ncm_annotation::set_vessel_sorting(ui.sortVessels->isChecked());

  // Zoom options
  cp.label_image_zoom_ = static_cast<ncm_qseries_preferences::ImageZoom>
      (ui.labelImageZoom->currentIndex());
  cp.add_vessels_zoom_ = static_cast<ncm_qseries_preferences::ImageZoom>
      (ui.addVesselsZoom->currentIndex());
  cp.set_vessel_props_zoom_ = static_cast<ncm_qseries_preferences::ImageZoom>
      (ui.setVesselPropsZoom->currentIndex());
  cp.add_haemorrhages_zoom_ = static_cast<ncm_qseries_preferences::ImageZoom>
      (ui.addHaemorrhagesZoom->currentIndex());

  // Apex labelling zoom factors
  cp.normal_vessel_zoom_ = ui.normalWindowSize->value();
  cp.enlarged_zoom_relative_ = ui.enlargedZoomRel->value();
  cp.giant_zoom_relative_ = ui.giantZoomRel->value();

  // Vessel drawing y-offset
  cp.vessel_path_yshift_ = ui.pathVerticalShift->value() * 0.01;

  // Grid size
  const double gridSpacingMm = ui.gridSpacingSpin->value();
  cp.set_grid_spacing_mm(gridSpacingMm);
  const double gridPixels = ui.gridPixelsEdit->text().toDouble();
  const double gridMillimetres = ui.gridMillimetresEdit->text().toDouble();
  double gridPixelsPerMm = gridPixels / gridMillimetres;
  cp.set_grid_pixels_per_mm(gridPixelsPerMm);

  // Warnings
  cp.warn_if_ungraded = ui.warnUngraded->isChecked();
  cp.warn_if_vessels_missing = ui.warnVesselsMissing->isChecked();
  cp.warn_if_sizes_missing = ui.warnSizesMissing->isChecked();
  v_app.warn_on_overwrite_size = ui.warnSizeOverwrite->isChecked();
  cp.warn_if_shapes_missing = ui.warnShapesMissing->isChecked();
  v_app.warn_on_overwrite_shape = ui.warnShapeOverwrite->isChecked();
  cp.warn_if_apices_missing = ui.warnApicesMissing->isChecked();
  cp.warn_if_paths_missing = ui.warnPathsMissing->isChecked();
  cp.warn_if_haemorrhages_missing = ui.warnHaemorrhagesMissing->isChecked();
}
