#ifndef DMK_LISTENER_H
#define DMK_LISTENER_H

// Listener.h: interface for the CListener class.
// 
// The CListener class is derived from GrabberListener. It overwrites
// the "frameReady()" method. In the frameReady method, the method
// "saveImage()" is called.
// "saveImage()" saves the specified buffer to a BMP file and calls a "Sleep(250)" 
// to simulate time consuming image processing. "saveImage()" is also called
// by the main() function of this example to save all buffers that have
// not been processed in the frameReady method.
//
// This class also overwrites the overlayCallback method to draw a 
// frame counter.
//
// The CListener object is registered with the parameter 
// eFRAMEREADY|eOVERLAYCALLBACK . 
//
//////////////////////////////////////////////////////////////////////

#include <QObject>
#include <QtCore>
#include <stdlib.h>
#include <stdio.h>
#include "tisudshl.h"

#define MESSAGEDEVICELOST WM_USER+90


using namespace DShowLib;

class CListener : public QObject, public GrabberListener   
{
	Q_OBJECT

	public:
		//void SetViewCWnd( CWnd *pView);
		CListener();
		virtual ~CListener();
		virtual void deviceLost( Grabber& param);
		virtual void frameReady( Grabber& param, smart_ptr<MemBuffer> mem_buffer, DWORD FrameNumber);

	signals:
		void frame_ready( smart_ptr<MemBuffer> mem_buffer, QTime curr_time);
		//void frame_ready();

	protected:
		
		//CWnd* m_pDrawCWnd;
		SIZE m_WindowSize;		// Size of the window in which to draw the buffer.
		void DrawBuffer( smart_ptr<MemBuffer> pBuffer);
		void DoImageProcessing( smart_ptr<MemBuffer> pBuffer);
		
				
};

#endif //DMK_LISTENER_H
