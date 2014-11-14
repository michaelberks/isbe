#ifndef NAILFOLD_QCAPTURE_SEARCH_H
#define NAILFOLD_QCAPTURE_SEARCH_H

#include <QDialog>
#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlQueryModel>
#include <QSortFilterProxyModel>
#include "ui_ncm_qcapture_search.h"
#include <nailfold/qtools/ncm_subject.h>

class ncm_qcapture_search : public QDialog
{
  Q_OBJECT

public:
	
	enum search_type {
		ImageSession,
		ReviewSession
	};

  ncm_qcapture_search(ncm_subject & subject, QString & study, const QSqlDatabase & database, const search_type search_type,
												QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture_search();



public slots:

protected:

private slots:

void on_chooseStudyComboBox_currentIndexChanged (const QString &text);
void on_chooseSubjectEdit_textChanged();
void on_selectSubjectPushButton_clicked();
void subject_selection_changed();
void closeEvent(QCloseEvent *ev);



private:
  //  Methods

	//Update GUI when text in choose subject box is changed
	void update_subject_filter();

	//Update display details of when a new subject is detected
	void update_subject_details();

  //: User interface object
  Ui::searchDialog ui;

	//Database object
	QSqlDatabase const &database_;

	//Currently selected study
	QString &selected_study_name_;
	int selected_study_id_;

	//Current subject name (in text edit box)
	QString subject_name_;

	//Current subject (loaded from SQL when subject selected)
	ncm_subject &subject_;

	//Type of session we're preforming a search for
	const search_type search_type_;

	//: Are we allowed to create new subjects? Yes for an imaging session,
	// no for reviewing/processing/reporting session
	const bool allow_new_subject_;

	//Model for subjects
	QSqlQueryModel *subjectsModel_;
	QSortFilterProxyModel *subjectFilterModel_;
};

#endif // NAILFOLD_QCAPTURE_SEARCH_H
