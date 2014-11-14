#include <nailfold/qtools/ncm_video_frame_queue.h>

ncm_video_frame_queue::ncm_video_frame_queue()
{
}

//: Add an item to the tail of the queue.
void ncm_video_frame_queue::enqueue(QSharedPointer<ncm_video_frame> qs)
{
	QMutexLocker ml(&queue_mutex_);
	//queue_mutex_.lock();
	  queue_.enqueue(qs);
	//queue_mutex_.unlock();
}

//: Remove and return the first item from the head of the queue.
QSharedPointer<ncm_video_frame> ncm_video_frame_queue::dequeue()
{
	QMutexLocker ml(&queue_mutex_);
	//queue_mutex_.lock();
	  QSharedPointer<ncm_video_frame> qs = queue_.dequeue();
	//queue_mutex_.unlock();

	return qs;
}

void ncm_video_frame_queue::empty()
{
	QMutexLocker ml(&queue_mutex_);
	queue_.clear();
}

//: Return a reference to the item at the head of the queue.
QSharedPointer<ncm_video_frame>& ncm_video_frame_queue::head()
{
	QMutexLocker ml(&queue_mutex_);
  return queue_.head();
}

//: Return a reference to the item at the tail of the queue.
QSharedPointer<ncm_video_frame>& ncm_video_frame_queue::tail()
{
	QMutexLocker ml(&queue_mutex_);
  return queue_.last();
}

//: Return true if queue is empty.
bool ncm_video_frame_queue::is_empty() const
{
  return queue_.isEmpty();
}
 
//: Return number of items in the queue.
int ncm_video_frame_queue::count() const
{
  return queue_.count();
}

