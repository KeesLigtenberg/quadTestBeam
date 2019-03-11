//rootmacro
//calculate chip position from Fred's measurements


/*
	 ______ ______
 	|      |      |
 	|      |      |
	|f____l|f____l|
	|		 GUARD		|
	|______ ______|
	|l    f|l    f|
 	|      |      |
	|______|______|
	

	Y	
	|__X

*/


//print vector
std::ostream& operator<<(std::ostream& out, const TVector2& v) {
	return out<<v.X()<<", "<<v.Y();
};


//calculate deviation from expected position in top left corner, taking into account possible rotations
void calculateDiffAndAngle(const TVector2& BLdiff, const TVector2& BRdiff) {
	auto BR=BRdiff+TVector2{255*55,0}; //BR coordinates
	auto BL=BLdiff; //origin is at BL
	auto BL_BR=BR-BL; //vector from BL->BR
	double angle=BL_BR.Phi();
	std::cout<<"origin offset by "<<BL<<" and angle is offset by "<< (angle > TMath::Pi() ? angle-2*TMath::Pi() : angle ) <<"\n";
	
}


void calculateAlignmentFromMeasurement() {
	TVector2 firstPad[4]={ {-3,20}, {18,3}, {11,21}, {5,8} };
	TVector2 lastPad[4]={ {2,7}, {22, 7}, {17,12}, {8,3} };
	
	for(int i=0; i<4; i++) {
		std::cout<<"Chip "<<i<<" ";
		calculateDiffAndAngle(firstPad[i], lastPad[i]);
	}

}
