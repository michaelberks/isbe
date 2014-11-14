#include "ncm_exposure_spinbox.h"
#include <vcl_iostream.h>
#include <vcl_cmath.h>

QExposureSpinBox::QExposureSpinBox( QWidget *parent )
		: QDoubleSpinBox( parent )
{
}

QExposureSpinBox::~QExposureSpinBox()
{
}

QString QExposureSpinBox::textFromValue ( double value ) const
{
	double new_val = A_ * exp(lambda_ * value);

	if ( new_val > 1 )
	{
		return QDoubleSpinBox::textFromValue( new_val );
	}
	else
	{
		QString str1 = "1 / ";
		QString str2;
		str2.setNum(1 / new_val, 'f', 2);
		str1.append( str2 );
		return str1;
	}

}

double QExposureSpinBox::valueFromText( const QString & text ) const
{
	vcl_cout << "Sub-class valueFromText called" << vcl_endl;

	double value = QDoubleSpinBox::valueFromText( text );
	return log(value / A_) / lambda_;
}

double QExposureSpinBox::valueFromExposure( double exposure )
{
	return log(exposure / A_) / lambda_;
}

double QExposureSpinBox::exposureFromValue( double value )
{
	return A_ * exp(lambda_ * value);
}

void QExposureSpinBox::setExposure( double value )
{
	QDoubleSpinBox::setValue( log(value / A_) / lambda_ );
}

void QExposureSpinBox::setExposureMin( double value )
{
	exposure_min_ = value;
	update_params();
}

void QExposureSpinBox::setExposureMax( double value )
{
	exposure_max_ = value;
	update_params();
}

void QExposureSpinBox::update_params()
{
	double x1 = this->minimum();
	double x2 = this->maximum();
	double y1 = exposure_min_;
	double y2 = exposure_max_;
	
	A_ = pow(y1, x2 / (x2-x1)) * 
       pow(y2, x1 / (x1-x2)); // = y1^(x2/(x2-x1)) * y2^(x1/(x1-x2));

	lambda_ = log(y1/y2) / (x1-x2);
}