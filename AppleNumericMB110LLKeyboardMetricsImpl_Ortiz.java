package keyboard;
/*Amanda Ortiz 
 * Spring 2016 - COSC 3327 - Algorithms - KART
 * Keyboard Metrics
 */

import static keyboard.Key.*;

import static keyboard.KeyLayout.COLEMAK;
import static keyboard.KeyLayout.DVORAK;
import static keyboard.KeyLayout.QWERTY;
import static keyboard.KeyLayout.ROTATION_13;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import combinatorics.Permutation;
import combinatorics.Permutation_Skeleton;

/**
 * @author skeleton
 *
 */
public class AppleNumericMB110LLKeyboardMetricsImpl_Ortiz implements KeyboardMetrics, AssignmentMetaData {
	private List<Key> vertexLabels;
	private int[][] adjacencyMatrix;
	private int[][] distanceMatrix;
	private Key homeKey;

	private static Map<KeyLayout, Key> keyLayoutToHomeKeyMap;
	private static Map<KeyLayout, Map<Key, Set<Key>>> keyLayoutToKeyToNeighborMapMap;

	static
	{
		keyLayoutToHomeKeyMap = new HashMap<KeyLayout, Key>();
		keyLayoutToHomeKeyMap.put(QWERTY, J);
		keyLayoutToHomeKeyMap.put(DVORAK, H);
		keyLayoutToHomeKeyMap.put(COLEMAK, N);
		keyLayoutToHomeKeyMap.put(ROTATION_13, W);

		keyLayoutToKeyToNeighborMapMap = new HashMap<KeyLayout, Map<Key, Set<Key>>>();
		Map<Key, Set<Key>> keyToNeighborMap_QWERTY = getKeyToNeighborMap_QWERTY();
		Map<Key, Set<Key>> keyToNeighborMap_DVORAK = getKeyToNeighborMap_Dvorak();
		Map<Key, Set<Key>> keyToNeighborMap_COLEMAK = getKeyToNeighborMap_Colemak();
		Map<Key, Set<Key>> keyToNeighborMap_ROT_13 = getKeyToNeighborMap_Rot13();

		keyLayoutToKeyToNeighborMapMap.put(QWERTY, keyToNeighborMap_QWERTY);
		keyLayoutToKeyToNeighborMapMap.put(DVORAK, keyToNeighborMap_DVORAK);
		keyLayoutToKeyToNeighborMapMap.put(COLEMAK, keyToNeighborMap_COLEMAK);
		keyLayoutToKeyToNeighborMapMap.put(ROTATION_13, keyToNeighborMap_ROT_13);
	}

	public AppleNumericMB110LLKeyboardMetricsImpl_Ortiz(KeyLayout keyLayout)
	{
		this.homeKey = keyLayoutToHomeKeyMap.get(keyLayout);
		Map<Key, Set<Key>> keyToNeighborsMap = keyLayoutToKeyToNeighborMapMap.get(keyLayout);
		init(keyToNeighborsMap, new ArrayList<Key>(keyToNeighborsMap.keySet()));
	}

	public void init(Map<Key, Set<Key>> physicalKeyToNeighborsMap, List<Key> vertexLabels)
	{
		this.vertexLabels = vertexLabels;
		this.adjacencyMatrix = getAdjacencyMatrix(physicalKeyToNeighborsMap, vertexLabels);
		this.distanceMatrix = getDistanceMatrix(adjacencyMatrix);
	}

	private static int[][] getAdjacencyMatrix(Map<Key, Set<Key>> physicalKeyToNeighborsMap, List<Key> vertexLabels)
	{
		assert physicalKeyToNeighborsMap.keySet().equals(new HashSet<Key>(vertexLabels)) : "vertexLabels inconsistent with physicalKeyToNeighborsMap! : vertexLabels = " + vertexLabels + " physicalKeyToNeighborsMap.keySet() = " + physicalKeyToNeighborsMap.keySet();
		final int SIZE = physicalKeyToNeighborsMap.keySet().size(); 
		int[][] adjacencyMatrix = new int[SIZE][SIZE];

		for(int row = 0; row < SIZE; row++)
		{
			for(int column = 0; column < SIZE; column++)
			{
				if(physicalKeyToNeighborsMap.get(vertexLabels.get(row)).contains(vertexLabels.get(column)))
				{
					adjacencyMatrix[row][column] = 1;
				}
				else
				{
					adjacencyMatrix[row][column] = 0;
				}
			}

		}



		return adjacencyMatrix;
	}

	//Matrix multiplication
	private static int[][] multiply(int[][] A, int[][] B)
	{
		int rowCount_A = A.length;
		assert rowCount_A > 0 : "rowCount_A = 0!";
		int columnCount_A = A[0].length;
		int rowCount_B = B.length;
		assert rowCount_B > 0 : "rowCount_B = 0!";
		int columnCount_B = B[0].length;
		assert columnCount_A == rowCount_B : "columnCount_A = " + columnCount_A + " <> " + rowCount_B + " = rowCount_B!";

		int[][] C = new int[rowCount_A][columnCount_B];
		for (int i = 0; i < rowCount_A; i++)
			for (int j = 0; j < columnCount_B; j++)
				for (int k = 0; k < columnCount_A; k++)
					C[i][j] += A[i][k] * B[k][j];

		return C;
	}

	private static int[][] getDistanceMatrix(int[][] adjacencyMatrix)
	{
		int vertexCount = adjacencyMatrix.length;
		assert vertexCount > 0 : "rowCount = 0!";
		boolean isMatrixFilled=false;
		int multiplicity = 1;
		int [] [] tempMatrix = adjacencyMatrix;
		int [] [] distanceMatrix = new int [vertexCount] [vertexCount];
		while (!isMatrixFilled) {
			for (int row=0; row<vertexCount; row++) {
				for (int column=0; column<vertexCount; column++) {
					if(tempMatrix[row][column]>0.0 && distanceMatrix[row][column]==0.0) {
						distanceMatrix[row][column]= multiplicity;
					}
				}
			}
			isMatrixFilled=true;
			for(int row=0; row<vertexCount; row++){
				for (int column=0; column<vertexCount; column++) {
					if (distanceMatrix[row][column]==0) {
						isMatrixFilled=false;
					}
				}
			}
			multiplicity++;
			tempMatrix=multiply(tempMatrix, adjacencyMatrix);

		}

		return distanceMatrix;

	}

	/* (non-Javadoc)
	 * @see keyboard.KeyboardMeasurements#getDistance(keyboard.PhysicalKey, keyboard.PhysicalKey)
	 */
	@Override
	public double getDistance(Key key1, Key key2) {
		int index1 = getIndex(vertexLabels, key1);
		int index2 = getIndex(vertexLabels, key2);
		if (key1==key2) {
			return 0;
		}
		return distanceMatrix[index1][index2];
	}

	private static <E> int getIndex(List<E> list, E element)
	{
		boolean foundIndex = false;
		int i = 0;
		while(!foundIndex && i < list.size())
		{
			foundIndex = (list.get(i) == element);
			if(!foundIndex) i++;
		}
		int rv = -1;
		if(foundIndex) rv = i;
		return rv;
	}

	@Override
	//******************//
	public double getDistance(String str) {

		double distance = 0;
		Key currentKey = homeKey;
		Set<Key> currentCharSet;
		Set <Key> nextCharSet;
		List <Key> currentCharList= new ArrayList<Key>();
		List <Key> nextCharList = new ArrayList<Key>();
		currentCharSet = getKeySet(str.charAt(0));
		currentCharList.addAll(currentCharSet);

		if (currentKey!=currentCharList.get(0)) {
			distance=getDistance(currentKey, currentCharList.get(0));
		}
		currentCharList.clear();

		for (int i=0; i<str.length()-1;i++) {
			currentCharSet = getKeySet(str.charAt(i));
			currentCharList.addAll(currentCharSet);
			nextCharSet = getKeySet(str.charAt(i+1));
			nextCharList.addAll(nextCharSet);

			if(currentCharList.get(0)==Key.M && str.charAt(i+1)==' ') {
				nextCharList.clear();
				nextCharList.add(Key.SPACEBAR_5);

			}
			if(str.charAt(i)==' ' && nextCharList.get(0)==Key.M) {

				currentCharList.clear();
				currentCharList.add(Key.SPACEBAR_5);
			}
			if(currentCharList.get(0)==Key.N && str.charAt(i+1)==' ') {
				nextCharList.clear();
				nextCharList.add(Key.SPACEBAR_4);
			}
			if(nextCharList.get(0)==Key.N  && str.charAt(i)==' ') {
				currentCharList.clear();
				currentCharList.add(Key.SPACEBAR_4);
			}
			if(currentCharList.get(0)==Key.B && str.charAt(i+1)==' ') {
				nextCharList.clear();
				nextCharList.add(Key.SPACEBAR_3);
			}
			if(nextCharList.get(0)==Key.B  && str.charAt(i)==' ') {
				currentCharList.clear();
				currentCharList.add(Key.SPACEBAR_3);
			}
			if(currentCharList.get(0)==Key.V && str.charAt(i+1)==' ') {
				nextCharList.clear();
				nextCharList.add(Key.SPACEBAR_2);
			}
			if(nextCharList.get(0)==Key.V  && str.charAt(i)==' ') {
				currentCharList.clear();
				currentCharList.add(Key.SPACEBAR_2);
			}			
			if(currentCharList.get(0)==Key.C && str.charAt(i+1)==' ') {
				nextCharList.clear();
				nextCharList.add(Key.SPACEBAR_1);
			}
			if(nextCharList.get(0)==Key.C  && str.charAt(i)==' ') {
				currentCharList.clear();
				currentCharList.add(Key.SPACEBAR_1);
			}

			distance+= getDistance(currentCharList.get(0), nextCharList.get(0));
			currentCharList.clear();
			nextCharList.clear();
		}	

		return distance;
	}

	private Key getClosestKey(Set<Key> keySet, Key key)
	{
		double minDistance = 0.0;
		List<Key> keyList = new ArrayList<Key>(keySet);
		Key minDistanceKey = null;

		//DO SOMETHING HERE...
		//getDistance() is involved...

		return minDistanceKey;
	}

	private static Set<Key> getKeySet(char character)
	{
		List<Key> keyList = Arrays.asList(Key.values());
		Set<Key> characterProducingKeysSet = new HashSet<Key>();
		for(int i = 0; i < keyList.size(); i++)
		{
			Key key = keyList.get(i);
			assert key != null : "key is null!";
			boolean keyProducesCharacter = (key.getNormalCharacter() != null && key.getNormalCharacter() == character) || (key.getShiftModifiedCharacter() != null && key.getShiftModifiedCharacter() == character);
			if(keyProducesCharacter) characterProducingKeysSet.add(key);
		}
		return characterProducingKeysSet;
	}

	private static Map<Key, Set<Key>> getKeyToNeighborMap_QWERTY()
	{
		Map<Key, Set<Key>> keyToNeighborSetMap = new HashMap<Key, Set<Key>>();
		keyToNeighborSetMap.put(Key.BACKTICK, getSet(TAB, ONE));
		keyToNeighborSetMap.put(Key.TAB, getSet(ONE, BACKTICK,Q));
		keyToNeighborSetMap.put(Key.ONE, getSet(TAB, BACKTICK, Q, TWO));
		keyToNeighborSetMap.put(Key.Q, getSet(TAB, ONE, TWO, W,A));
		keyToNeighborSetMap.put(Key.SHIFT_1, getSet(A, Z));
		keyToNeighborSetMap.put(Key.A, getSet(SHIFT_1, Z, Q, W, S));
		keyToNeighborSetMap.put(Key.Z, getSet(SHIFT_1, A, X, S));
		keyToNeighborSetMap.put(Key.TWO, getSet(ONE, Q, W, THREE));
		keyToNeighborSetMap.put(Key.W, getSet(TWO, Q, A, S, E, THREE));
		keyToNeighborSetMap.put(Key.S, getSet(W, A, Z, X, D, E));
		keyToNeighborSetMap.put(Key.X, getSet(Z, S, D, C));
		keyToNeighborSetMap.put(Key.THREE, getSet(TWO, W, E, FOUR));
		keyToNeighborSetMap.put(Key.E, getSet(THREE, W, S, D, R, FOUR));
		keyToNeighborSetMap.put(Key.D, getSet(E, R, F, C, X, S));
		keyToNeighborSetMap.put(Key.C, getSet(X, D, F, V, SPACEBAR_1));
		keyToNeighborSetMap.put(Key.SPACEBAR_1, getSet(Key.C));
		keyToNeighborSetMap.put(Key.FOUR, getSet(THREE, E, R, FIVE));
		keyToNeighborSetMap.put(Key.R, getSet(FOUR, E, D, F, T, FIVE));
		keyToNeighborSetMap.put(Key.F, getSet(R, D, C, V, G, T));
		keyToNeighborSetMap.put(Key.V, getSet(C, F, G, B, SPACEBAR_2));
		keyToNeighborSetMap.put(Key.SPACEBAR_2, getSet(Key.V));
		keyToNeighborSetMap.put(Key.FIVE, getSet(FOUR, R, T, SIX));
		keyToNeighborSetMap.put(Key.T, getSet(FIVE, R, F, G, Y, SIX));
		keyToNeighborSetMap.put(Key.G, getSet(T, F, V, B, H, Y));
		keyToNeighborSetMap.put(Key.B, getSet(V, G, H, N, SPACEBAR_3));
		keyToNeighborSetMap.put(Key.SPACEBAR_3, getSet(B));
		keyToNeighborSetMap.put(Key.SIX, getSet(FIVE, T, Y, SEVEN));
		keyToNeighborSetMap.put(Key.Y, getSet(SIX, T, G, H, U, SEVEN));
		keyToNeighborSetMap.put(Key.H, getSet(Y, G, B, N, J, U));
		keyToNeighborSetMap.put(Key.N, getSet(B, H, J, M, SPACEBAR_4));
		keyToNeighborSetMap.put(Key.SPACEBAR_4, getSet(Key.N));
		keyToNeighborSetMap.put(Key.SEVEN, getSet(SIX, Y, U, EIGHT));
		keyToNeighborSetMap.put(Key.U, getSet(SEVEN, Y, H, J, I, EIGHT));
		keyToNeighborSetMap.put(Key.J, getSet(U, H, N, M, K, I));
		keyToNeighborSetMap.put(Key.M, getSet(N, J, K, COMMA, SPACEBAR_5));
		keyToNeighborSetMap.put(Key.SPACEBAR_5, getSet(Key.M));
		keyToNeighborSetMap.put(Key.EIGHT, getSet(SEVEN, U, I, NINE));
		keyToNeighborSetMap.put(Key.I, getSet(EIGHT, U, J, K, O, NINE));
		keyToNeighborSetMap.put(Key.K, getSet(I, J, M, COMMA, L, O));
		keyToNeighborSetMap.put(Key.COMMA, getSet(M, K, L, PERIOD));
		keyToNeighborSetMap.put(Key.NINE, getSet(EIGHT, I, O, ZERO));
		keyToNeighborSetMap.put(Key.O, getSet(NINE, I, K, L, P, ZERO));
		keyToNeighborSetMap.put(Key.L, getSet(O, K, COMMA, PERIOD, SEMICOLON, P));
		keyToNeighborSetMap.put(Key.PERIOD, getSet(COMMA, L, SEMICOLON, FORESLASH));
		keyToNeighborSetMap.put(Key.ZERO, getSet(NINE, O, P, MINUS));
		keyToNeighborSetMap.put(Key.P, getSet(ZERO, O, L, SEMICOLON, LEFT_BRACKET, MINUS));
		keyToNeighborSetMap.put(Key.SEMICOLON, getSet(P, L, PERIOD, FORESLASH, TICK, LEFT_BRACKET));
		keyToNeighborSetMap.put(Key.FORESLASH, getSet(PERIOD, SEMICOLON, TICK, SHIFT_2));
		keyToNeighborSetMap.put(Key.MINUS, getSet(ZERO, P, LEFT_BRACKET, EQUALS));
		keyToNeighborSetMap.put(Key.LEFT_BRACKET, getSet(MINUS, P, SEMICOLON, TICK, RIGHT_BRACKET, EQUALS));
		keyToNeighborSetMap.put(Key.TICK, getSet(LEFT_BRACKET, SEMICOLON, FORESLASH, SHIFT_2, RETURN, RIGHT_BRACKET));
		keyToNeighborSetMap.put(Key.SHIFT_2, getSet(FORESLASH, TICK, RETURN));
		keyToNeighborSetMap.put(Key.EQUALS, getSet(MINUS, LEFT_BRACKET, RIGHT_BRACKET));
		keyToNeighborSetMap.put(Key.RIGHT_BRACKET, getSet(EQUALS, LEFT_BRACKET, TICK, RETURN, BACKSLASH));
		keyToNeighborSetMap.put(Key.RETURN, getSet(BACKSLASH, RIGHT_BRACKET, TICK, SHIFT_2));
		keyToNeighborSetMap.put(Key.BACKSLASH, getSet(RIGHT_BRACKET, RETURN));


		return keyToNeighborSetMap;
	}

	private static Map<Key, Set<Key>> getKeyToNeighborMap_Dvorak()
	{
		Map<Key, Set<Key>> keyToNeighborSetMap = new HashMap<Key, Set<Key>>();
		keyToNeighborSetMap.put(Key.SPACEBAR_1, getSet(Key.J));
		keyToNeighborSetMap.put(Key.SPACEBAR_2, getSet(Key.K));
		keyToNeighborSetMap.put(Key.SPACEBAR_3, getSet(B));
		keyToNeighborSetMap.put(Key.SPACEBAR_4, getSet(Key.B));
		keyToNeighborSetMap.put(Key.SPACEBAR_5, getSet(Key.M));
		//bottom row
		keyToNeighborSetMap.put(Key.SHIFT_1, getSet(SEMICOLON, A));
		keyToNeighborSetMap.put(Key.SEMICOLON, getSet(Q, O, A, SHIFT_1));
		keyToNeighborSetMap.put(Key.Q, getSet(J, E, O, SEMICOLON));
		keyToNeighborSetMap.put(Key.J, getSet(K, U, E, Q, SPACEBAR_1));
		keyToNeighborSetMap.put(Key.K, getSet(X, I, U, J, SPACEBAR_2));
		keyToNeighborSetMap.put(Key.X, getSet(B, D, I, K, SPACEBAR_3));
		keyToNeighborSetMap.put(Key.B, getSet(M, H, D, X, SPACEBAR_4));
		keyToNeighborSetMap.put(Key.M, getSet(W, T, H, B, SPACEBAR_5));
		keyToNeighborSetMap.put(Key.W, getSet(V, N, T, M));
		keyToNeighborSetMap.put(Key.V, getSet(Z, S, N, W));
		keyToNeighborSetMap.put(Key.Z, getSet(SHIFT_2, MINUS, S, V));
		keyToNeighborSetMap.put(Key.SHIFT_2, getSet(RETURN, MINUS, Z));
		//second row
		keyToNeighborSetMap.put(Key.A, getSet(O, COMMA, TICK, SHIFT_1, SEMICOLON));
		keyToNeighborSetMap.put(Key.O, getSet(E, PERIOD, COMMA, A, SEMICOLON, Q));
		keyToNeighborSetMap.put(Key.E, getSet(U, P, PERIOD, O, Q, J));
		keyToNeighborSetMap.put(Key.U, getSet(I, Y, P, E, J, K));
		keyToNeighborSetMap.put(Key.I, getSet(D, F, Y, U, K, X));
		keyToNeighborSetMap.put(Key.D, getSet(H, G, F, I,X,B));
		keyToNeighborSetMap.put(Key.H, getSet(T, C, G, D,B,M));
		keyToNeighborSetMap.put(Key.T, getSet(N, R, C, H, M, W));
		keyToNeighborSetMap.put(Key.N, getSet(S, L, R, T, W, V));
		keyToNeighborSetMap.put(Key.S, getSet(MINUS, FORESLASH, L, N, V, Z));
		keyToNeighborSetMap.put(Key.MINUS, getSet(RETURN, EQUALS, FORESLASH, S, Z, SHIFT_2));
		keyToNeighborSetMap.put(Key.RETURN, getSet(BACKSLASH, EQUALS, MINUS, SHIFT_2));
		//third row
		keyToNeighborSetMap.put(Key.TAB, getSet(TICK, ONE, BACKTICK));
		keyToNeighborSetMap.put(Key.TICK, getSet(COMMA, TWO, ONE, TAB, A));
		keyToNeighborSetMap.put(Key.COMMA, getSet(PERIOD, THREE, TWO, TICK, A, O));
		keyToNeighborSetMap.put(Key.PERIOD, getSet(P, FOUR, THREE, COMMA, P, E));
		keyToNeighborSetMap.put(Key.P, getSet(Y, FIVE, FOUR, PERIOD, E, U));
		keyToNeighborSetMap.put(Key.Y, getSet(F, SIX, FIVE, P, U, I));
		keyToNeighborSetMap.put(Key.F, getSet(G, SEVEN, SIX, Y, I, D));
		keyToNeighborSetMap.put(Key.G, getSet(C, EIGHT, SEVEN, F, D, H));
		keyToNeighborSetMap.put(Key.C, getSet(R, NINE, EIGHT, G, H, T));
		keyToNeighborSetMap.put(Key.R, getSet(L, ZERO, NINE, C, T, N));
		keyToNeighborSetMap.put(Key.L, getSet(FORESLASH, LEFT_BRACKET, ZERO, R, N, S));
		keyToNeighborSetMap.put(Key.FORESLASH, getSet(EQUALS, RIGHT_BRACKET, LEFT_BRACKET, L, S, MINUS));
		keyToNeighborSetMap.put(Key.EQUALS, getSet(BACKSLASH, RIGHT_BRACKET, FORESLASH, MINUS, RETURN));
		keyToNeighborSetMap.put(Key.BACKSLASH, getSet(EQUALS, RETURN));
		//fourth row
		keyToNeighborSetMap.put(Key.BACKTICK, getSet(ONE, TAB));
		keyToNeighborSetMap.put(Key.ONE, getSet(TWO, BACKTICK, TAB, TICK));
		keyToNeighborSetMap.put(Key.TWO, getSet(THREE, ONE, TICK, COMMA));
		keyToNeighborSetMap.put(Key.THREE, getSet(FOUR, TWO, COMMA, PERIOD));
		keyToNeighborSetMap.put(Key.FOUR, getSet(FIVE, THREE, PERIOD, P));
		keyToNeighborSetMap.put(Key.FIVE, getSet(SIX, FOUR, P, Y));
		keyToNeighborSetMap.put(Key.SIX, getSet(SEVEN, FIVE, Y, F));
		keyToNeighborSetMap.put(Key.SEVEN, getSet(EIGHT, SIX, F, G));
		keyToNeighborSetMap.put(Key.EIGHT, getSet(NINE, SEVEN, G, C));
		keyToNeighborSetMap.put(Key.NINE, getSet(ZERO, EIGHT, C, R));
		keyToNeighborSetMap.put(Key.ZERO, getSet(LEFT_BRACKET, NINE, R, L));
		keyToNeighborSetMap.put(Key.LEFT_BRACKET, getSet(RIGHT_BRACKET, ZERO, L, FORESLASH));
		keyToNeighborSetMap.put(Key.RIGHT_BRACKET, getSet(LEFT_BRACKET, FORESLASH, EQUALS));
		return keyToNeighborSetMap;
	}

	private static Map<Key, Set<Key>> getKeyToNeighborMap_Colemak()
	{
		Map<Key, Set<Key>> keyToNeighborSetMap = new HashMap<Key, Set<Key>>();

		//spacebars
		keyToNeighborSetMap.put(Key.SPACEBAR_1, getSet(C));
		keyToNeighborSetMap.put(Key.SPACEBAR_2, getSet(V));
		keyToNeighborSetMap.put(Key.SPACEBAR_3, getSet(B));
		keyToNeighborSetMap.put(Key.SPACEBAR_4, getSet(K));
		keyToNeighborSetMap.put(Key.SPACEBAR_5, getSet(M));
		//bottom row
		keyToNeighborSetMap.put(Key.SHIFT_1, getSet(Z, A));
		keyToNeighborSetMap.put(Key.Z, getSet(X, R, A, SHIFT_1));
		keyToNeighborSetMap.put(Key.X, getSet(C, S, R, Z));
		keyToNeighborSetMap.put(Key.C, getSet(V, T, S, X, SPACEBAR_1));
		keyToNeighborSetMap.put(Key.V, getSet(B, D, T, C, SPACEBAR_2));
		keyToNeighborSetMap.put(Key.B, getSet(K, H, D, V, SPACEBAR_3));
		keyToNeighborSetMap.put(Key.K, getSet(M, N, H, B, SPACEBAR_4));
		keyToNeighborSetMap.put(Key.M, getSet(COMMA, E, N, K, SPACEBAR_5));
		keyToNeighborSetMap.put(Key.COMMA, getSet(PERIOD, I, E, M));
		keyToNeighborSetMap.put(Key.PERIOD, getSet(FORESLASH, O, I, COMMA));
		keyToNeighborSetMap.put(Key.FORESLASH, getSet(SHIFT_2, TICK, O, PERIOD));
		keyToNeighborSetMap.put(Key.SHIFT_2, getSet(RETURN, TICK, FORESLASH));

		//second row
		keyToNeighborSetMap.put(Key.A, getSet(R, W, Q, SHIFT_1, Z));
		keyToNeighborSetMap.put(Key.R, getSet(S, F, W, A, Z, X));
		keyToNeighborSetMap.put(Key.S, getSet(T, P, F, R, X, C));
		keyToNeighborSetMap.put(Key.T, getSet(D, G, P, S, C, V));
		keyToNeighborSetMap.put(Key.D, getSet(H, J, G, T, V, B));
		keyToNeighborSetMap.put(Key.H, getSet(N, L, J, D, B, K));
		keyToNeighborSetMap.put(Key.N, getSet(E, U, L, H, K, M));
		keyToNeighborSetMap.put(Key.E, getSet(I, Y, U, N, M, COMMA));
		keyToNeighborSetMap.put(Key.I, getSet(O, SEMICOLON, Y, E, COMMA, PERIOD));
		keyToNeighborSetMap.put(Key.O, getSet(TICK, LEFT_BRACKET, SEMICOLON, I, PERIOD, FORESLASH, SHIFT_2));
		keyToNeighborSetMap.put(Key.TICK, getSet(RETURN, RIGHT_BRACKET, LEFT_BRACKET, O, FORESLASH));
		keyToNeighborSetMap.put(Key.RETURN, getSet(BACKSLASH, RIGHT_BRACKET, TICK, SHIFT_2));
		//third row
		keyToNeighborSetMap.put(Key.TAB, getSet(Q, ONE, BACKTICK));
		keyToNeighborSetMap.put(Key.Q, getSet(W, TWO, ONE, TAB, A));
		keyToNeighborSetMap.put(Key.W, getSet(F, THREE, TWO, Q, A, R));
		keyToNeighborSetMap.put(Key.F, getSet(P, FOUR, THREE, W, R, S));
		keyToNeighborSetMap.put(Key.P, getSet(G, FIVE, FOUR, F, S, T));
		keyToNeighborSetMap.put(Key.G, getSet(J, SIX, FIVE, P, T, D));
		keyToNeighborSetMap.put(Key.J, getSet(L, SEVEN, SIX, G, D, H));
		keyToNeighborSetMap.put(Key.L, getSet(U, EIGHT, SEVEN, J, H, N));
		keyToNeighborSetMap.put(Key.U, getSet(Y, NINE, EIGHT, L, N, E));
		keyToNeighborSetMap.put(Key.Y, getSet(SEMICOLON, ZERO, NINE, U, E, I));
		keyToNeighborSetMap.put(Key.SEMICOLON, getSet(LEFT_BRACKET, MINUS, ZERO, Y, I, O));
		keyToNeighborSetMap.put(Key.LEFT_BRACKET, getSet(RIGHT_BRACKET, EQUALS, MINUS, SEMICOLON, O, TICK));
		keyToNeighborSetMap.put(Key.RIGHT_BRACKET, getSet(BACKSLASH, EQUALS, LEFT_BRACKET, TICK, RETURN));
		keyToNeighborSetMap.put(Key.BACKSLASH, getSet(RIGHT_BRACKET, RETURN));
		//fourth row
		keyToNeighborSetMap.put(Key.BACKTICK, getSet(ONE, TAB));
		keyToNeighborSetMap.put(Key.ONE, getSet(TWO, Q, TAB, BACKTICK));
		keyToNeighborSetMap.put(Key.TWO, getSet(THREE, ONE, W, Q));
		keyToNeighborSetMap.put(Key.THREE, getSet(FOUR, TWO, COMMA, PERIOD));
		keyToNeighborSetMap.put(Key.FOUR, getSet(FIVE, THREE, P, F));
		keyToNeighborSetMap.put(Key.FIVE, getSet(SIX, FOUR, G, P));
		keyToNeighborSetMap.put(Key.SIX, getSet(SEVEN, FIVE, J, G));
		keyToNeighborSetMap.put(Key.SEVEN, getSet(EIGHT, SIX, L, J));
		keyToNeighborSetMap.put(Key.EIGHT, getSet(NINE, SEVEN, U, L));
		keyToNeighborSetMap.put(Key.NINE, getSet(ZERO, EIGHT, Y, U));
		keyToNeighborSetMap.put(Key.ZERO, getSet(MINUS, NINE, SEMICOLON, Y));
		keyToNeighborSetMap.put(Key.MINUS, getSet(EQUALS, ZERO, LEFT_BRACKET, SEMICOLON));
		keyToNeighborSetMap.put(Key.EQUALS, getSet(MINUS, BACKSLASH, RIGHT_BRACKET));

		return keyToNeighborSetMap;
	}

	private static Map<Key, Set<Key>> getKeyToNeighborMap_Rot13()
	{
		Map<Key, Set<Key>> keyToNeighborSetMap = new HashMap<Key, Set<Key>>();
		keyToNeighborSetMap.put(Key.BACKTICK, getSet(TAB, ONE));
		keyToNeighborSetMap.put(Key.TAB, getSet(ONE, BACKTICK,D));
		keyToNeighborSetMap.put(Key.ONE, getSet(TAB, BACKTICK, D, TWO));
		keyToNeighborSetMap.put(Key.D, getSet(TAB, ONE, TWO, N, J));
		keyToNeighborSetMap.put(Key.SHIFT_1, getSet(N, M));
		keyToNeighborSetMap.put(Key.N, getSet(SHIFT_1, M, F, J, D));
		keyToNeighborSetMap.put(Key.M, getSet(SHIFT_1, N, K, F));
		keyToNeighborSetMap.put(Key.TWO, getSet(ONE, D, J, THREE));
		keyToNeighborSetMap.put(Key.J, getSet(TWO, D, N, F, R, THREE));
		keyToNeighborSetMap.put(Key.F, getSet(N, M, K, Q, R, J));
		keyToNeighborSetMap.put(Key.K, getSet(M, F, Q, P));
		keyToNeighborSetMap.put(Key.THREE, getSet(TWO, J, R, FOUR));
		keyToNeighborSetMap.put(Key.R, getSet(THREE, J, F, Q, E, FOUR));
		keyToNeighborSetMap.put(Key.Q, getSet(R, F, K, P, S, E));
		keyToNeighborSetMap.put(Key.P, getSet(K, Q, S, I, SPACEBAR_1));
		keyToNeighborSetMap.put(Key.SPACEBAR_1, getSet(Key.P));
		keyToNeighborSetMap.put(Key.FOUR, getSet(THREE, E, R, FIVE));
		keyToNeighborSetMap.put(Key.E, getSet(FOUR, R, Q, S, G, FIVE));
		keyToNeighborSetMap.put(Key.S, getSet(Q, P, I, T, G, E));
		keyToNeighborSetMap.put(Key.I, getSet(P, S, T, O, SPACEBAR_2));
		keyToNeighborSetMap.put(Key.SPACEBAR_2, getSet(Key.I));
		keyToNeighborSetMap.put(Key.FIVE, getSet(FOUR, E, G, SIX));
		keyToNeighborSetMap.put(Key.G, getSet(FIVE, E, S, T, L, SIX));
		keyToNeighborSetMap.put(Key.T, getSet(G, S, I, O, U, L));
		keyToNeighborSetMap.put(Key.O, getSet(I, T, U, A, SPACEBAR_3));
		keyToNeighborSetMap.put(Key.SPACEBAR_3, getSet(O));
		keyToNeighborSetMap.put(Key.SIX, getSet(FIVE, G, L, SEVEN));
		keyToNeighborSetMap.put(Key.L, getSet(SIX, G, T, U, H, SEVEN));
		keyToNeighborSetMap.put(Key.U, getSet(L, T, O, A, W, H));
		keyToNeighborSetMap.put(Key.A, getSet(O,U, W, Z, SPACEBAR_4));
		keyToNeighborSetMap.put(Key.SPACEBAR_4, getSet(Key.A));
		keyToNeighborSetMap.put(Key.SEVEN, getSet(SIX, L, H, EIGHT));
		keyToNeighborSetMap.put(Key.H, getSet(SEVEN, L, U, W, V, EIGHT));
		keyToNeighborSetMap.put(Key.W, getSet(H, U, A, Z, X, V));
		keyToNeighborSetMap.put(Key.Z, getSet(A, W, X, COMMA, SPACEBAR_5));
		keyToNeighborSetMap.put(Key.SPACEBAR_5, getSet(Key.Z));
		keyToNeighborSetMap.put(Key.EIGHT, getSet(SEVEN, H, V, NINE));
		keyToNeighborSetMap.put(Key.V, getSet(EIGHT, H, W, X, B, NINE));
		keyToNeighborSetMap.put(Key.X, getSet(V, W, Z, COMMA, Y, B));
		keyToNeighborSetMap.put(Key.COMMA, getSet(Z, X, Y, PERIOD));
		keyToNeighborSetMap.put(Key.NINE, getSet(EIGHT, V, B, ZERO));
		keyToNeighborSetMap.put(Key.B, getSet(NINE, V, X, Y, C, ZERO));
		keyToNeighborSetMap.put(Key.Y, getSet(C, B, COMMA, PERIOD, SEMICOLON, X));
		keyToNeighborSetMap.put(Key.PERIOD, getSet(COMMA, Y, SEMICOLON, FORESLASH));
		keyToNeighborSetMap.put(Key.ZERO, getSet(NINE, B, V, MINUS));
		keyToNeighborSetMap.put(Key.C, getSet(ZERO, B, Y, SEMICOLON, LEFT_BRACKET, MINUS));
		keyToNeighborSetMap.put(Key.SEMICOLON, getSet(Y, C, PERIOD, FORESLASH, TICK, LEFT_BRACKET));
		keyToNeighborSetMap.put(Key.FORESLASH, getSet(PERIOD, SEMICOLON, TICK, SHIFT_2));
		keyToNeighborSetMap.put(Key.MINUS, getSet(ZERO, C, LEFT_BRACKET, EQUALS));
		keyToNeighborSetMap.put(Key.LEFT_BRACKET, getSet(MINUS, C, SEMICOLON, TICK, RIGHT_BRACKET, EQUALS));
		keyToNeighborSetMap.put(Key.TICK, getSet(LEFT_BRACKET, SEMICOLON, FORESLASH, SHIFT_2, RETURN, RIGHT_BRACKET));
		keyToNeighborSetMap.put(Key.SHIFT_2, getSet(FORESLASH, TICK, RETURN));
		keyToNeighborSetMap.put(Key.EQUALS, getSet(MINUS, LEFT_BRACKET, RIGHT_BRACKET));
		keyToNeighborSetMap.put(Key.RIGHT_BRACKET, getSet(EQUALS, LEFT_BRACKET, TICK, RETURN, BACKSLASH));
		keyToNeighborSetMap.put(Key.RETURN, getSet(BACKSLASH, RIGHT_BRACKET, TICK, SHIFT_2));
		keyToNeighborSetMap.put(Key.BACKSLASH, getSet(RIGHT_BRACKET, RETURN));


		return keyToNeighborSetMap;
	}



	private static Set<Key> getSet(Key... keys)
	{
		return new HashSet<Key>(Arrays.asList(keys));
	}

	@Override
	public String getFirstNameOfSubmitter() {
		// TODO Auto-generated method stub
		return "Amanda";
	}

	@Override
	public String getLastNameOfSubmitter() {
		
		return "Ortiz";
	}

	@Override
	public double getHoursSpentWorkingOnThisAssignment() {
		
		return 35;
	}

	@Override
	public int getScoreAgainstTestCasesSubset() {
		
		return 100;
	}


}