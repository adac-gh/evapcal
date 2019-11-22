	PROGRAM CASCADE
C++
C	Main program calling CASCIP, which contains the original program
C	CASCADE by Puhlhofer (further comments in CASCIP)
C	Two versions of the subroutine CASCIP are available, with file names
C	CASCIPM.FOR and CASCNM.FOR. The first one ne.

	COMMON/CASCMIN/ INDCASC
	INDCASC=1
	CALL CASCIP

	STOP
	END
