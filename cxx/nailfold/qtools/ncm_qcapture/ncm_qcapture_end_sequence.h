#ifndef NAILFOLD_qcapture_end_sequence_H
#define NAILFOLD_qcapture_end_sequence_H

#include <QDialog>

#include "ui_ncm_qcapture_end_sequence.h"
#include <nailfold/qtools/ncm_image_sequence.h>

class ncm_qcapture_end_sequence : public QDialog
{
  Q_OBJECT

public:
  ncm_qcapture_end_sequence(ncm_image_sequence & sequence, const QString study_name, const QString subject_study_name,
												QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture_end_sequence();

public slots:

protected:

private slots:

	void set_digit(int index);
	void on_sequenceNotesPlainTextEdit_textChanged();
	void on_qualityCheckBox_toggled(bool checked);
	void on_sequenceNameLineEdit_textChanged();
	void on_endSequenceButton_clicked();
	//void on_


private:
  //  Methods
	void createHandRadioGroup();

  //: User interface object
  Ui::end_sequenceDialog ui;

	//: Radio group to select dominant hand
	QButtonGroup *handRadioGroup_;

	//Current sequence
	ncm_image_sequence & sequence_;
};

#endif // NAILFOLD_qcapture_end_sequence_H
