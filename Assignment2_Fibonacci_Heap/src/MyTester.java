import java.util.Arrays;
import java.util.*;
public class MyTester {
    public static void main(String[] args) {
        FibonacciHeap heap2 = new FibonacciHeap();
        int p = 25;
        int a = (int)Math.pow(2,p);
        int m = (int)(Math.log(a) / Math.log(2));
        long start1 = System.currentTimeMillis();
        FibonacciHeap.HeapNode[] nodeArray = new FibonacciHeap.HeapNode[a+1];
        for (int i = a -1; i>= -1; i--){
            nodeArray[i+1] = heap2.insert(i);        }
        heap2.deleteMin();
        for (int j = m; j>0; j--){
            heap2.decreaseKey(nodeArray[a - (int)Math.pow(2,j) + 2], a + 1);
        }
        long start2 = System.currentTimeMillis();
        System.out.println(FibonacciHeap.totalLinks());
        System.out.println(FibonacciHeap.totalCuts());
        System.out.println(heap2.potential());
        System.out.println(start2-start1);

        /*test1();
        testSameNodes();*/
    }
    public static void testSameNodes() {
        FibonacciHeap heap = new FibonacciHeap();
        List<FibonacciHeap.HeapNode> lst = new ArrayList<>();
        samenodesInsert(heap,lst);
        nodesInsert(heap,lst);
        int size=lst.size();
        for (int i = 0; i <lst.size()-1 ; i++) {
            deletex(heap,lst);
            heap.checkParent();
            if(heap.findMin()==null)
            {
                System.out.println("errorrrrrrr");
            }
            check(heap);

        }


    }
    public static void test1() {
        FibonacciHeap heap = new FibonacciHeap();
        List<FibonacciHeap.HeapNode> lst = new ArrayList<>();
        nodesInsert(heap,lst);
        deletex(heap,lst);
        decreasekeys(heap,lst);
        System.out.println(Arrays.toString(heap.countersRep()));
        check(heap);

    }
    public static void test() {
        FibonacciHeap heap = new FibonacciHeap();
        List<Integer> lst = new ArrayList<>();
        heapInsert(heap, lst);
        int size = lst.size();
        for (int i = 0; i < size; i++) {
            boolean a = check(heap);
            if (a == false) {
                return;
            }

            delete(heap, lst);

        }
    }

    public static void heapInsert(FibonacciHeap heap, List<Integer> arr) {
        for (int i = 0; i < 1000; i++) {
            heap.insert(i);
            arr.add(i);
        }

    }

    public static void delete(FibonacciHeap heap, List<Integer> arr) {
        Random index = new Random();
        arr.remove(0);
        heap.deleteMin();

    }

    public static boolean check(FibonacciHeap heap) {
        boolean a = heap.checkRank();
        if (a == false) {
            System.out.println("check rank is not good");
            return false;
        }
        a = heap.checkCounters();
        if (a == false) {
            System.out.println("counters rep not good");
            return false;
        }
        a = heap.checkParent();
        if (a == false) {
            System.out.println("parent is not good");
            return false;
        }

        return true;
    }

    public static void nodesInsert(FibonacciHeap heap, List<FibonacciHeap.HeapNode> arr) {
        for (int i = 0; i < 1000; i++) {
            FibonacciHeap.HeapNode a = new FibonacciHeap.HeapNode(i);
            heap.insertNode(a);
            arr.add(a);
        }

    }
    public static void samenodesInsert(FibonacciHeap heap, List<FibonacciHeap.HeapNode> arr) {
        for (int i = 0; i <1000 ; i++) {
            FibonacciHeap.HeapNode a = new FibonacciHeap.HeapNode(1);
            heap.insertNode(a);
            arr.add(a);
        }
    }

    public static void decreasekeys(FibonacciHeap heap, List<FibonacciHeap.HeapNode> arr) {
        for (int i = 0; i < arr.size(); i++) {
            FibonacciHeap.HeapNode curr = arr.get(i);
            if (curr.getParent() != null) {
                heap.decreaseKey(curr, (curr.getKey() - curr.getParent().getKey() + 10));

            }

        }
    }
        public static void deletex(FibonacciHeap heap,List<FibonacciHeap.HeapNode> arr)
        {
            FibonacciHeap.HeapNode a=arr.get(0);
            arr.remove(0);
            heap.delete(a);

        }

}






