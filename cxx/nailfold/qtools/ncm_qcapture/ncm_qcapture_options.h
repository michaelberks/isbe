#ifndef NAILFOLD_QCAPTURE_OPTIONS_H
#define NAILFOLD_QCAPTURE_OPTIONS_H

#include <QDialog>
#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlQueryModel>
#include <QSortFilterProxyModel>
#include "ui_ncm_qcapture_options.h"
#include <nailfold/qtools/ncm_subject.h>
#include "ncm_qcapture_preferences.h"

class ncm_qcapture_options : public QDialog
{
  Q_OBJECT

public:
  ncm_qcapture_options(ncm_qcapture_preferences & preferences,  
		QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture_options();

	//Define an enum if changed
	enum {
		PreferencesChanged = 0,
		PreferencesUnchanged = 1};

public slots:

protected:

private slots:
	//3 main controls
	void on_okButton_clicked();
	void on_cancelButton_clicked();
	void on_applyButton_clicked();

	//---------------------------------------------------------
	//General tab
	void on_autoLoadCamera_toggled(bool checked);
	void on_rotateFrames_toggled(bool checked);
	void on_deleteEmptySessions_toggled(bool checked);
	void on_processUnacceptable_toggled(bool checked);

	void on_saveDirSelect_clicked();
	void on_sequenceDataLineEdit_textChanged(QString text);
	void on_framePrefixLineEdit_textChanged(QString text);

	void on_maxSaveFramesSpinBox_valueChanged(int value);
	void on_maxMosaicFramesSpinBox_valueChanged(int value);
	void on_maxMosaicSizeSpinBox_valueChanged(int value);
	void on_nthAlignFrameSpinBox_valueChanged(int value);

	//----------------------------------------------------------
	//Motors tab
	void on_homeMotorsOnStartup_toggled(bool checked);

	//Set speeds of motors for X/Y panning
	void on_dblclickSpeedMax_valueChanged(double value);
	void on_joystickSpeedMax_valueChanged(double value);

	//Set speed/dist for Z motor
	void on_zarrowDist_valueChanged(double value);
	void on_zscrollDist_valueChanged(double value);
	
	//Calibrate the pixel resolution
	void on_pixelResImageSpinBox_valueChanged(int value);
	void on_pixelResMotorSpinBox_valueChanged(int value);

	//----------------------------------------------------------
	//Auto focus tab
	void on_computeSharpness_toggled(bool checked);
	void on_enableAutofocus_toggled(bool checked);
	void on_nSharpnessSpinBox_valueChanged(int value);
	void on_modulateWeight_valueChanged(int value);

	

//void on_chooseStudyComboBox_currentIndexChanged (const QString &text);
//void on_chooseSubjectEdit_textChanged();
//void on_selectSubjectPushButton_clicked();
//void subject_selection_changed();



private:
  //  Methods
	void write_values();
	void preferences_modified();

	// Variables

	Ui::optionsDialog ui;

	//New preferences option modified during dialog
	ncm_qcapture_preferences new_preferences_;

	//Current preferences on opening dialog - won't be modified unless 'apply' or 'ok' selected
	ncm_qcapture_preferences & current_preferences_;

	int preferences_changed_;
};

#endif // NAILFOLD_QCAPTURE_OPTIONS_H
