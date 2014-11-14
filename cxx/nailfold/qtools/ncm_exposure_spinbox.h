#ifndef NCM_QCAPTURE_EXPOSURE_SPINBOX
#define NCM_QCAPTURE_EXPOSURE_SPINBOX

#include <QtCore>
#include <QDoubleSpinBox>
#include <QListView>
#include <vcl_iostream.h>

class QExposureSpinBox : public QDoubleSpinBox
{
    Q_OBJECT

	public:
		QExposureSpinBox( QWidget *parent = 0 );
		~QExposureSpinBox();
		
		virtual QString	textFromValue ( double value ) const;
		virtual double valueFromText( const QString & text ) const;
		double valueFromExposure( double exposure );
		double exposureFromValue( double value );

		void setExposure( double value );
		void setExposureMin( double value );
		void setExposureMax( double value );

	private:
		//Exposure min and max
		double exposure_min_;
		double exposure_max_;
		double A_;
		double lambda_;
		void update_params();
};

#endif //NCM_QCAPTURE_EXPOSURE_SPINBOX