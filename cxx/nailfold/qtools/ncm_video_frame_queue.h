#ifndef NCM_VIDEO_FRAME_QUEUE_H
#define NCM_VIDEO_FRAME_QUEUE_H

#include <QtCore>
#include "ncm_video_frame.h"

class ncm_video_frame_queue
{
// INTERFACE

public:

	ncm_video_frame_queue();

	void enqueue(QSharedPointer<ncm_video_frame> qs);
	QSharedPointer<ncm_video_frame> dequeue();
	void empty();

  //: Return a reference to the item at the head of the queue.
	QSharedPointer<ncm_video_frame>& head();

  //: Return a reference to the item at the tail of the queue.
	QSharedPointer<ncm_video_frame>& tail();

  //: Return true if queue is empty.
  bool is_empty() const;

  //: Return number of items in queue.
  int count() const;


// IMPLEMENTATION

protected:

private: // Variables
	QQueue< QSharedPointer<ncm_video_frame> > queue_;
	QMutex queue_mutex_;
};

#endif // NCM_VIDEO_FRAME_QUEUE_H