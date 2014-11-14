#ifndef __ncm_hand_h__
#define __ncm_hand_h__

//:
// \file
// \brief Define enums of hand (left/right and digits)
// \author Mike Berks
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_string.h>

class ncm_hand
{
public:
	enum Hand { NotKnown = 0,
							Left = 1,
							Right = 2};

	enum Digit { None = 0,
							 D1 = 1,
							 D2 = 2,
							 D3 = 3,
							 D4 = 4,
							 D5 = 5};

	static Hand fromString(vcl_string str)
	{
		if(str == "left") return Left; 

		if(str == "Left") return Left;

		if(str == "l") return Left;

		if(str == "L") return Left;	

		if(str == "right") return Right; 

		if(str == "Right") return Right;

		if(str == "r") return Right;

		if(str == "R") return Right;
			
		//Otherwise we don't know what to do...
		return NotKnown;
	}

		static vcl_string toString(Hand hand)
		{
			switch(hand)
			{
				case (Left):
					return "Left";
				case (Right):
					return "Right";
				case (NotKnown):
					return "Not known";
				default:
					return "Not known";
			}
		}

		static vcl_string toLetter(Hand hand)
		{
			switch(hand)
			{
				case (Left):
					return "L";
				case (Right):
					return "R";
				case (NotKnown):
					return "N";
				default:
					return "N";
			}
		}

		static vcl_string toLetter(Digit digit)
		{
			switch(digit)
			{
				case (D1):
					return "1";
				case (D2):
					return "2";
				case (D3):
					return "3";
				case (D4):
					return "4";
				case (D5):
					return "5";
				case (None):
					return "N";
				default:
					return "N";
			}
		}

};

#endif // __ncm_hand_h__