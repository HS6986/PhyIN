/* This is a snippet of the Java code for PhyIN v. 1.0 as implemented in the pre-4.0 development version of Mesquite (www.mesquiteproject.org).

This snippet resides in the class mesquite.molec.FlagByPhyIN.FlagByPhyIN 

It is available in the development branch of Mesquite (https://github.com/MesquiteProject/MesquiteCore/tree/development) */


	double proportionIncompat = proportionIncompatDEFAULT; //(-p)
	int blockSize = blockSizeDEFAULT; //(-b)
	int neighbourDistance = neighbourDistanceDEFAULT; //(-d)
	MesquiteBoolean treatGapAsState = new MesquiteBoolean(treatGapAsStateDEFAULT); //(-e)
	MesquiteBoolean examineOnlySelectedTaxa = new MesquiteBoolean(examineOnlySelectedTaxaDEFAULT); //(no equivalent in python script)
	
	boolean[] hasConflict;
	boolean[] toSelect;
	int[] taxonSequenceStart, taxonSequenceEnd;
	int NUMSTATES = 4;
	int numTaxaWithSequence = 0;
	boolean[][] statePairs;

	/*....................................................................................*/
	int getEffectiveState(long stateSet) {
		if (treatGapAsState.getValue() && CategoricalState.isInapplicable(stateSet))
			return NUMSTATES-1; 
		if (CategoricalState.isCombinable(stateSet) && !CategoricalState.hasMultipleStates(stateSet))
			return CategoricalState.getOnlyElement(stateSet);
		return -1;
	}
	int getFirstInSequence(CategoricalData data, int it) {
		for (int ic=0; ic<data.getNumChars(); ic++)
			if (!data.isInapplicable(ic, it))
				return ic;
		return -1;
	}
	int getLastInSequence(CategoricalData data, int it) {
		for (int ic=data.getNumChars()-1; ic>=0; ic--)
			if (!data.isInapplicable(ic, it))
				return ic;
		return -1;
	}

	boolean anyOtherConnectForward(int i, int k) {
		for (int m =0; m<NUMSTATES;m++) 
			if (k != m && statePairs[i][m]) {
				return true;
			}
		return false;
	}
	boolean anyOtherConnectBackward(int i, int k) {
		for (int m =0; m<NUMSTATES;m++) 
			if (i != m && statePairs[m][k]) {
				return true;
			}
		return false;
	}
	boolean anyOfGraphLeft() {
		for (int i =0; i<NUMSTATES; i++) 
			for (int k =0; k<NUMSTATES;k++) {
				if (statePairs[i][k]) {
					return true;
				}
			}
		return false;
	}
	/*....................................................................................*/
	boolean trimSingleConnects() {
		//turn a statePairs[i][k] to false if it is the only pairing/edge among statePairs[i]
		for (int i =0; i<NUMSTATES; i++) {
			for (int k =0; k<NUMSTATES;k++) {
				if (statePairs[i][k]) { //edge is in graph
					if (!anyOtherConnectForward(i, k)) {
						statePairs[i][k] = false;
						return true;
					}
					if (!anyOtherConnectBackward(i, k)) {
						statePairs[i][k] = false;
						return true;
					}
				}
			}
		}
		return false;
	}
	/*....................................................................................*/
	boolean anyTaxaSelected = false;
	boolean includeTaxon(int it, CategoricalData data) {
		return !anyTaxaSelected || (!examineOnlySelectedTaxa.getValue() || data.getTaxa().getSelected(it));
	}
	/*....................................................................................*/
	boolean areIncompatible(CategoricalData data, int ic, int ic2) {
		if (ic>=data.getNumChars() || ic <0 || ic2>=data.getNumChars() || ic2 <0)
			return false;

		for (int i =0; i<NUMSTATES; i++) for (int k =0; k<NUMSTATES;k++) statePairs[i][k]= false;
	
		//first, harvest all patterns between the two columns
		for (int it = 0; it < data.getNumTaxa(); it++) {
			// only look at taxa for which ic and ic2 are within their sequence (i.e. not in terminal gap region)
			if (includeTaxon(it, data) && taxonSequenceStart[it]>=0 && ic>= taxonSequenceStart[it] && ic <= taxonSequenceEnd[it] && ic2>= taxonSequenceStart[it] && ic2 <= taxonSequenceEnd[it]) {
				int state1 = getEffectiveState(data.getState(ic, it));
				int state2 = getEffectiveState(data.getState(ic2, it));
				if (state1 >=0 && state1<NUMSTATES && state2 >=0 && state2<NUMSTATES) {
					statePairs[state1][state2] = true;
				}
			}
		}

		//Test of compatibility: Look for cycles in state to state occupancy graph (M. Steel)
		int stopper = 1000; //merely to prevent infinite loop in case of bug
		while (trimSingleConnects() && (stopper--) > 0) {
			if (stopper == 1)
				MesquiteMessage.warnUser("ALERT SOMETHING WENT WRONG WITH THE PhyIN calculations (trimSingleConnects too many iterations)");
		}
		// if anything is left of the graph, then it's incompatible
		if (anyOfGraphLeft()) 
			return true;

		return false;
	}

	/*....................................................................................*/
	/** Marks those characters in conflict with others*/
	public void markThoseWithConflict(CategoricalData data){
		if (taxonSequenceStart == null || taxonSequenceStart.length != data.getNumTaxa()) {
			taxonSequenceStart = new int[data.getNumTaxa()];
			taxonSequenceEnd = new int[data.getNumTaxa()];
		}
		anyTaxaSelected = data.getTaxa().anySelected();
		if (hasConflict == null || hasConflict.length != data.getNumChars()) {
			hasConflict = new boolean[data.getNumChars()];
			toSelect = new boolean[data.getNumChars()];
		}
		for (int ic=0; ic<hasConflict.length; ic++) {
			hasConflict[ic] = false;
			toSelect[ic] = false;
		}
		numTaxaWithSequence = 0;
		for (int it=0; it<data.getNumTaxa(); it++) {
			if (includeTaxon(it, data)){
				taxonSequenceStart[it] = getFirstInSequence(data, it);
				taxonSequenceEnd[it] = getLastInSequence(data, it);
				if (taxonSequenceStart[it]>=0)
					numTaxaWithSequence++;
			}
		}

		NUMSTATES = data.getMaxState()+1; 
		if (treatGapAsState.getValue())
			NUMSTATES++;//one extra for gaps as state
		if (statePairs == null || statePairs.length != NUMSTATES)
			statePairs = new boolean[NUMSTATES][NUMSTATES];

		for (int ic=0; ic<data.getNumChars(); ic++) {
			for (int d = 1; d<neighbourDistance+1; d++)
				if (areIncompatible(data, ic, ic+d)) {
					hasConflict[ic] = true;
					if (ic+d<hasConflict.length)
						hasConflict[ic+d] = true;
				}
		}
	}
	/*....................................................................................*/
	void selectSpanByProportion(int ic, boolean[] hasConflict, boolean[] toSelect, double proportionIncomp, int spanSize) {
		if (!hasConflict[ic])
			return;
		int count = 0;
		int lastHit = -1;
		for (int k = 0; k<spanSize && ic+k<hasConflict.length; k++) {
			if (hasConflict[ic+k]) {
				count++;
				lastHit = ic+k;
			}
		}
		if (1.0*count/spanSize >= proportionIncomp) {
			for (int k = ic; k<=lastHit && k<toSelect.length; k++) {
				toSelect[k]=true;
			}
		}

	}

	/*....................................................................................*/
	public MatrixFlags flagMatrix(CharacterData data, MatrixFlags flags) {
		if (data!=null && data.getNumChars()>0 && data instanceof CategoricalData){
			if (flags == null)
				flags = new MatrixFlags(data);
			else {
				flags.reset(data);
			}
			markThoseWithConflict((CategoricalData)data);
			for (int ic=0; ic<hasConflict.length; ic++) {
				selectSpanByProportion(ic, hasConflict, toSelect, proportionIncompat, blockSize);
			}

			for (int ic=0; ic<data.getNumChars(); ic++) {
				if (toSelect[ic])
					flags.setCharacterFlag(ic, true);
			}
		}
		return flags;
	}