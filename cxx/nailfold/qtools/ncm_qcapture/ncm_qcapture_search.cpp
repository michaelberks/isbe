#include "ncm_qcapture_search.h"

// VXL libraries
#include <vcl_iostream.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>
#include <nailfold/qtools/ncm_qvesselitem.h>
#include <nailfold/ncm_hand.h>
#include "ncm_qcapture_new_subject.h"

// Qt libraries
#include <QMessageBox>
#include <QColorDialog>
#include <QtSql/QSqlRecord>
//
//: Constructor
ncm_qcapture_search::ncm_qcapture_search(ncm_subject &subject, QString &study, const QSqlDatabase & database, const search_type search_type,
																					QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags),
	subject_(subject),
	selected_study_name_(study),
	database_(database),
  subject_name_(""),
	search_type_(search_type),
	allow_new_subject_(search_type==ImageSession)
{
  // setup the UI
  ui.setupUi(this);

	//Initialise subjects view model
	subjectsModel_ = new QSqlQueryModel;
	subjectFilterModel_ = new QSortFilterProxyModel();

	//Set the selection mode of the subject view
	ui.existingSujectsView->setSelectionMode(QAbstractItemView::SingleSelection);

	//Set the subject details box to be read only
	ui.selectedSubjectDetails->setReadOnly(true);

	
	QSqlQuery query(database_);

	//Get studies stored in database and populate combobox
	ui.chooseStudyComboBox->clear();

	query.exec("SELECT study_name FROM study");
  while (query.next()) {
      QString study = query.value(0).toString();
			ui.chooseStudyComboBox->addItem(study);
	}
	selected_study_name_ = ui.chooseStudyComboBox->currentText();

	//Get users stored in database and populate combobox
	ui.imagedByComboBox->clear();

	query.exec("SELECT username FROM user");
  while (query.next()) {
    QString user = query.value(0).toString();
		ui.imagedByComboBox->addItem(user);
  }

}

//
//: Destructor
ncm_qcapture_search::~ncm_qcapture_search()
{
}

//
//  Public slots
//

//
//: What to do when the user hits 'OK'

//
//  Private slots
//
void ncm_qcapture_search::on_chooseStudyComboBox_currentIndexChanged (const QString &text)
{

	//Store the selected study
	selected_study_name_ = text;

	QString study_query_str = "SELECT ncm_study_id FROM study WHERE study_name = \'";
	study_query_str += selected_study_name_;
	study_query_str += "\'";
	//vcl_cout << "Query string = " << study_query_str.toStdString() << vcl_endl;

	QSqlQuery query(database_);
	query.exec(study_query_str);
	query.first();
	selected_study_id_ = query.value(0).toInt();

	//vcl_cout << "Selected study " << selected_study_name_.toStdString() << vcl_endl;
	//vcl_cout << "Selected study ID " << selected_study_id_ << vcl_endl;

	//Get list of subjects from this study
	QString subject_query_str = "SELECT * FROM subject WHERE subject_study_id = ";
	subject_query_str += QString::number(selected_study_id_);

	//Query the database and get the column index of the subject name
	subjectsModel_->setQuery(subject_query_str, database_);
	int fieldNo = subjectsModel_->record().indexOf("subject_study_name");
	
	//Link the SQL model to the filter model so we can filter on subject name
	subjectFilterModel_->setSourceModel( subjectsModel_ );
	subjectFilterModel_->setFilterKeyColumn(fieldNo);

	//Link the filter model to the list view to display the filtered subjects
	ui.existingSujectsView->setModel(subjectFilterModel_);
	ui.existingSujectsView->setModelColumn(fieldNo);
  ui.existingSujectsView->show();

	//Connect up the subject view selection model so we process selection changes
	connect(ui.existingSujectsView->selectionModel(), SIGNAL( selectionChanged (QItemSelection, QItemSelection)),
				this, SLOT(subject_selection_changed()) );

	update_subject_filter();
}

//Update list of subjects as text is entered
void ncm_qcapture_search::on_chooseSubjectEdit_textChanged() 
{
	subject_name_ = ui.chooseSubjectEdit->text();
	update_subject_filter();
	

}

void ncm_qcapture_search::on_selectSubjectPushButton_clicked()
{
	if (search_type_ == ReviewSession && subject_.last_imaged().isNull())
	{
		QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("Subject has no image sessions to review. \n Please select another subject");
		msgBox.exec();
		return;
	}
	
	else if (allow_new_subject_ && ui.selectSubjectPushButton->text() == "Create new subject" )
	{
		//vcl_cout << "Create new subject" << vcl_endl;
		
		//Open a new subject dialog window - pass in the subject object to have its fields completed by the user
		ncm_qcapture_new_subject new_subject_ui(subject_, selected_study_id_, selected_study_name_, subject_name_, this);
		new_subject_ui.exec();

		QSqlQuery query(database_);
		bool success = subject_.sql_write(query); //What to do if we can't write to the database?
	}

	//else - Nothing to do in this case, the subject is as we want and already in the database

	done(0); //We're done! close dialog with return code 0 (success)
	
}

void ncm_qcapture_search::subject_selection_changed()
{
	//Is a valid subject selected?
	if ( ui.existingSujectsView->selectionModel()->hasSelection() )
	{
		//Get the subject record from the sql model, then load this record into the current subject object
		QModelIndex filterIndex = ui.existingSujectsView->selectionModel()->selectedIndexes().first();
		QModelIndex srcIndex = subjectFilterModel_->mapToSource(filterIndex);
		QSqlQuery query(database_);
		subject_.sql_read( subjectsModel_->record(srcIndex.row()), query );

		//Update the information
		update_subject_details();

		//Enable the select button
		ui.selectSubjectPushButton->setEnabled(true);
	}
	else if( subjectFilterModel_->rowCount() > 0 )
	{
		//Selection has been filtered out, but we still have matching items - disable select button, clear current subject
		ui.selectSubjectPushButton->setEnabled(false);
		subject_.reset();
		ui.selectedSubjectDetails->setText("");
	}

	//else - No matching subjects, this case already dealt with in update_filter_selection

}

void ncm_qcapture_search::closeEvent(QCloseEvent *ev)
{
	subject_.reset();
}

//
//  Private functions
//
void ncm_qcapture_search::update_subject_filter()
{
	//Filter the subjects in the list view
	QString filter_str = "^" + subject_name_;// + "*";
	subjectFilterModel_->setFilterRegExp(filter_str);

	//What to do if no matching subjects
	if (subjectFilterModel_->rowCount() == 0 )
	{
		//Clear the current subject and the subject text box
		subject_.reset();
		ui.selectedSubjectDetails->setText("");

		if (allow_new_subject_)
		{
			//Update the behaviour of the "select case" button
			ui.selectSubjectPushButton->setText("Create new subject");
			ui.selectSubjectPushButton->setToolTip("Create a new subject for this imaging session");

			//if subject isn't empty, enable the button	
			ui.selectSubjectPushButton->setEnabled(!subject_name_.toStdString().empty());
		}
	}
	//We have subjects, again - change back to select
	else if (ui.selectSubjectPushButton->text() == "Create new subject" )
	{
		ui.selectSubjectPushButton->setText("Select subject");
		ui.selectSubjectPushButton->setToolTip("Choose selected subject for a new imaging session");
		ui.selectSubjectPushButton->setEnabled(false); //Disabled until a selection is made
	} 

}

void ncm_qcapture_search::update_subject_details()
{
	QString display_text = 
		"Study name: " + subject_.name_in_study() + "\n"
		"Subject ID: " + QString::number(subject_.ncm_id()) + "\n"
		"Dominant hand : " + QString(ncm_hand::toString(subject_.dominant_hand()).c_str()) + "\n"
		"Last imaged: ";

	QDateTime date = subject_.last_imaged();
	
	if (date.isNull())
		display_text += "(no previous sessions) \n";
	else
		display_text += (date.toString("ddd MMMM d yyyy") + "\n"); 

	display_text += ("Subject notes: " + subject_.notes());

	ui.selectedSubjectDetails->setText(display_text);
}