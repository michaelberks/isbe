#ifndef NCM_DMK_GRAB_FRAME_FILTER_H
#define NCM_DMK_GRAB_FRAME_FILTER_H

/**
 *	This frame filter applies a binarization on the image data.
 *	If enabled, every gray value greater of equal to a specified threshold is changed to 
 *	the maximum gray value, every other gray value is changed to zero.
 *
 *	Allowed input types: Y800, RGB8
 *
 *	Output types: Y800, RGB8, the input type determines the output type.
 *
 *	Parameters:
 *		enable:	Use enable( bool bEnable ) to enable or disable binarization.
 *				If binarization is disabled, the image data is not modified.
 *		threshold:
 *				Use setThreshold( int th ) to set the threshold for the
 *				binarization.
 *
 */
#include "tisudshl.h"
#include <QtCore>
#include <QSharedPointer.h>
#include <QObject>
#include <QTime>
#include <vil/vil_image_view.h>
#include <vcl_cstring.h>
#include <nailfold/qtools/ncm_qcapture/ncm_qcapture_data_manager.h>
#include <nailfold/qtools/ncm_video_frame.h>

class ncm_dmk_grab_frame_filter :	public QObject, public DShowLib::FrameFilterImpl<ncm_dmk_grab_frame_filter> 
{
	Q_OBJECT

public:
	ncm_dmk_grab_frame_filter();

	// FrameFilterImpl implementation
	static	DShowLib::FilterInfo getStaticFilterInfo();

	// IFrameFilter implementation
	virtual	void getSupportedInputTypes( DShowLib::FrameTypeInfoArray& arr ) const;
	virtual	bool getTransformOutputTypes( const DShowLib::FrameTypeInfo& in_type,
										  DShowLib::FrameTypeInfoArray& out_types ) const;
	virtual	bool transform( const DShowLib::IFrame& src, DShowLib::IFrame& dest );

	void set_queue(ncm_video_frame_queue* _queue);

signals:
		void frame_ready();

private:
	ncm_video_frame_queue* queue_;
};


#endif // NCM_DMK_GRAB_FRAME_FILTER_H