import java.util.*;

public class RandomTester {

    public static final Random RANDOM = new Random(1);
    public static final int NUMBER_OF_OPERATIONS = (int) Math.pow(2, 10);
    public static final int MAX_BOUND_FOR_K = 3000; //do max k inserts random
    public static final int DELTA_FOR_DECREASE_KEY = 1000; //we can replace it with random value or something
    public static final boolean VERBOSE = true;

    public static final int[] OPERATION_ON_EMPTY_HEAP = new int[]{0};
    public static final int[] OPERATION_ON_HEAP = new int[]{0,1,2,3,4,5,6}; //not using 7 for now
    public static final List<Integer> DE_LAZY_HEAP_OPERATION = Arrays.asList(2,3);

    // Operation map
    // 0 - insert Random
    // 1 - insert new minimum
    // 2 - delete random
    // 3 - delete minimum
    // 4 - insert k Random
    // 5 - meld
    // 6 - decreaseKey random - TODO check node after update
    // 7 - decreaseKey minimum - TODO check node after update

    // @pre keyToNodeMap.isEmpty() == false
    public static int getMinimumKeyInMap(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap) {
        int min = Integer.MAX_VALUE;
        for (int key : keyToNodeMap.keySet()) {
            if (key < min) {
                min = key;
            }
        }
        return min;
    }

    public static <T> void printErrorMessage(String test, T expectValue, T value) {
        if (VERBOSE) {
            System.out.printf("!%s, expected: %s, got: %s\n", test, expectValue, value);
        }
    }

    public static boolean isValidSize(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
        boolean isCorrectSize = keyToNodeMap.size() == fibonacciHeap.size();
        boolean isCorrectEmpty = keyToNodeMap.isEmpty() == fibonacciHeap.isEmpty();
        if (!isCorrectSize) { printErrorMessage("isCorrectSize", keyToNodeMap.size(), fibonacciHeap.size()); }
        if (!isCorrectEmpty) { printErrorMessage("isCorrectEmpty", keyToNodeMap.isEmpty(), fibonacciHeap.isEmpty()); }
        return isCorrectSize && isCorrectEmpty;
    }

    public static boolean isValidMin(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
        FibonacciHeap.HeapNode expectedMinNode = keyToNodeMap.isEmpty() ? null : keyToNodeMap.get(getMinimumKeyInMap(keyToNodeMap));
        FibonacciHeap.HeapNode minNode = fibonacciHeap.findMin();
        boolean isCorrectMinNode = expectedMinNode == minNode;
        if (!isCorrectMinNode) { printErrorMessage("isCorrectMinNode", expectedMinNode, minNode); }
        return isCorrectMinNode;
    }

    public static boolean isValidCountRep(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
        int size = keyToNodeMap.size();
        int[] expectedRep;
        fibonacciHeap.checkRank();
        if (size == 0) {
            expectedRep = new int[]{};
        } else {
            String binary = Integer.toBinaryString(size);
            expectedRep = new int[binary.length()];
            for (int index = 0; index < binary.length(); index++) {
                expectedRep[binary.length() - 1 - index] = binary.charAt(index) == '1' ? 1 : 0;
            }
        }
        int[] countersRep = fibonacciHeap.countersRep();
        boolean isCorrectCountersRep = Arrays.equals(expectedRep, countersRep);
        if (!isCorrectCountersRep) {
            printErrorMessage("isCorrectCountersRep", Arrays.toString(expectedRep), Arrays.toString(countersRep));}
        return true;
    }

    public static boolean isValidHeap(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap, boolean checkCountRep) {
        return isValidSize(keyToNodeMap, fibonacciHeap)
                && isValidMin(keyToNodeMap, fibonacciHeap)&& (!checkCountRep || isValidCountRep(keyToNodeMap, fibonacciHeap));
    }

    public static void test() {
        FibonacciHeap fibonacciHeap = new FibonacciHeap();
        HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap = new HashMap<>();
        int maxSize = 0;
        int countMeldOps = 0;
        int decKeyOps = 0;
        for (int steps = 0; steps < NUMBER_OF_OPERATIONS; steps++) {
            int operation = Utils.pickRandomOperation(keyToNodeMap);
            if (operation == 5) countMeldOps++;
            if (operation == 6 || operation == 7) decKeyOps++;
            if (VERBOSE) { Utils.printStatus(keyToNodeMap, fibonacciHeap); }
            Utils.performOperation(operation, keyToNodeMap, fibonacciHeap);
            boolean isValid = isValidHeap(keyToNodeMap, fibonacciHeap, DE_LAZY_HEAP_OPERATION.contains(operation));
            if (keyToNodeMap.size() > maxSize) maxSize = keyToNodeMap.size();
            if(!isValid) {
                System.out.println("not valid heap, failed after " + steps + " steps");
                break;
            }
        }
        System.out.println("max size the heap was at is: " + maxSize);
        System.out.println("total number of melds is: " + countMeldOps);
        System.out.println("total number of dec keys is: " + decKeyOps);
    }

    public static void main(String[] args) {
        test();
    }

    static class Utils {

        public static void insertNumber(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap, int key) {
            keyToNodeMap.put(key, fibonacciHeap.insert(key));
        }

        //I PUT NOT VERBOSE HERE TO AVOID INSANE SPAM
        public static void insertRandomNumber(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            while (true) {
                int key = RANDOM.nextInt();
                if (!keyToNodeMap.containsKey(key)) {
                    if (!VERBOSE) { System.out.println("insertRandomNumber key: " + key); }
                    Utils.insertNumber(keyToNodeMap, fibonacciHeap, key);
                    return;
                }
            }
        }

        public static void insertKRandomNumbers(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            for (int index = 0; index < MAX_BOUND_FOR_K; index++) {
                insertRandomNumber(keyToNodeMap, fibonacciHeap);
            }
        }

        public static void insertNewMinimum(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            int newMinimum = getMinimumKeyInMap(keyToNodeMap) - 1;
            if (VERBOSE) { System.out.println("insertNewMinimum key: " + newMinimum); }
            insertNumber(keyToNodeMap, fibonacciHeap, newMinimum);
        }

        public static int pickRandomKey(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap) {
            List<Integer> nodes = new ArrayList<>(keyToNodeMap.keySet());
            Collections.shuffle(nodes, RANDOM);
            int index = RANDOM.nextInt(nodes.size());
            return nodes.get(index);
        }

        public static void deleteKey(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap, int key) {
            FibonacciHeap.HeapNode node = keyToNodeMap.remove(key);
            fibonacciHeap.delete(node);
        }

        // @pre keyToNodeMap.isEmpty() == false
        public static void deleteRandomNode(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            int key = pickRandomKey(keyToNodeMap);
            if (VERBOSE) { System.out.println("deleteRandomNode key: " + key); }
            deleteKey(keyToNodeMap, fibonacciHeap, key);
        }

        // @pre keyToNodeMap.isEmpty() == false
        public static void deleteMinimumKey(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            int minimumKey = getMinimumKeyInMap(keyToNodeMap);
            if (VERBOSE) { System.out.println("deleteMinimumKey key: " + minimumKey); }
            deleteKey(keyToNodeMap, fibonacciHeap, minimumKey);
        }

        public static void meldRandomHeap(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            FibonacciHeap fibonacciHeap2 = new FibonacciHeap();
            insertKRandomNumbers(keyToNodeMap, fibonacciHeap2);
            if (VERBOSE) { System.out.printf("Melding random heap with %d elements\n", MAX_BOUND_FOR_K);}
            fibonacciHeap.meld(fibonacciHeap2);
        }

        public static void decreaseKey(int key, String msg, HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            FibonacciHeap.HeapNode node = keyToNodeMap.get(key);
            if (VERBOSE) {System.out.println(msg + key);}
            fibonacciHeap.decreaseKey(node, DELTA_FOR_DECREASE_KEY);
            keyToNodeMap.remove(key);
            keyToNodeMap.put(key - DELTA_FOR_DECREASE_KEY, node);
        }

        public static void decreaseRandomKey(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            int key = pickRandomKey(keyToNodeMap);
            decreaseKey(key, "decreaseRandomKey key: ", keyToNodeMap, fibonacciHeap);
        }

        public static void decreaseMinimumKey(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            int key = getMinimumKeyInMap(keyToNodeMap);
            decreaseKey(key, "decreaseMinimumKey key: ", keyToNodeMap, fibonacciHeap);
        }

        public static int pickRandomOperation(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap) {
            int upperBound;
            if (keyToNodeMap.isEmpty()) {
                upperBound = OPERATION_ON_EMPTY_HEAP.length;
            } else {
                upperBound = OPERATION_ON_HEAP.length;
            }
            return RANDOM.nextInt(upperBound);
        }

        public static void performOperation(int operation, HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            switch (operation) {
                case 0: {
                    insertRandomNumber(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 1: {
                    insertNewMinimum(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 2: {
                    deleteRandomNode(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 3: {
                    deleteMinimumKey(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 4: {
                    insertKRandomNumbers(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 5: {
                    meldRandomHeap(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 6: {
                    decreaseRandomKey(keyToNodeMap, fibonacciHeap);
                    break;
                }
                case 7: { //not using this for now (bugs with minValue)
                    decreaseMinimumKey(keyToNodeMap, fibonacciHeap);
                }
                default: {};
            }
        }

        public static void printStatus(HashMap<Integer, FibonacciHeap.HeapNode> keyToNodeMap, FibonacciHeap fibonacciHeap) {
            int size = keyToNodeMap.size();
            int minNodeKey = getMinimumKeyInMap(keyToNodeMap);
            int[] countRep = fibonacciHeap.countersRep();
            System.out.printf("size: %s, minNodeKey: %s, countRep: %s\n", size, minNodeKey, Arrays.toString(countRep) );
        }
    }
}