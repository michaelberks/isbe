#ifndef NAILFOLD_QFINAL_OPTIONS_H
#define NAILFOLD_QFINAL_OPTIONS_H

#include <QDialog>
#include <QDoubleValidator>

#include "ui_ncm_qfinal_options.h"

#include "ncm_qfinal_preferences.h"

class ncm_qfinal_options : public QDialog
{
  Q_OBJECT

public:
  ncm_qfinal_options(ncm_qfinal_preferences& preferences,
                      QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qfinal_options();

public slots:
  virtual void accept();

protected:

private slots:
  void on_apexColorButton_clicked();
  void on_vessel_colour_clicked();
  void on_vessel_colour_selected_clicked();
  void on_gridPixelsEdit_textEdited(const QString& text);
  void on_gridMillimetresEdit_textEdited(const QString& text);

private:
  //  Methods

  //: Get value from environment (i.e. get existing settings)
  void read_values();

  //: Copy values to environment (i.e. save new settings)
  void write_values();

  //: Helper functions for colour pushbuttons
  void apply_colour_to(QPushButton* const button, QRgb const& rgb);
  QRgb get_colour_from(QPushButton const* const button);

  //: Update labels of grid pixel->mm conversions
  void updateGridLabels();

  //  Variables

  //: Reference to the existing preferences. This is updated in write_values().
  ncm_qfinal_preferences& current_preferences_;

  //: Copy of the existing preferences. These are modified while the dialog is
  //  open, but not stored permanently until write_values().
  ncm_qfinal_preferences new_preferences_;

  //: User interface object
  Ui::optionsDialog ui;

  //: Map from rezoom factor labels to integer values
  QMap<QString, int> rezoom_factor_map_;

  QDoubleValidator doubleValidator_;

  QButtonGroup* userLevelGroup_;
};

#endif // NAILFOLD_QFINAL_OPTIONS_H
