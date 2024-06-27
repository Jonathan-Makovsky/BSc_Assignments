import java.util.ArrayList;// delete before submission
import java.util.Arrays;
import java.util.List;

/**
 * FibonacciHeap
 *
 * An implementation of a Fibonacci Heap over integers.
 */
public class FibonacciHeap {
    private HeapNode min;
    private HeapNode first;
    private HeapNode last;
    private int size;
    private static int totalCuts=0;
    private static int totalLinks=0;
    private int trees;
    private int marked;
    public FibonacciHeap() {
        this.min = null;
        this.first=null;
        this.last=null;
        this.size=0;
        this.trees=0;
        this.marked=0;
    }

    public FibonacciHeap(HeapNode node) {
        this.min = node;
        node.setKey(node.getKey());
        node.setRank(node.getRank());
        node.setMark(node.getMark());
        node.setChild(null);
        node.setNext(null);
        node.setPrev(null);
        node.setParent(null);
        this.first=node;
        this.last=node;
        this.size=1;
        this.trees=1;
        this.marked=0;
    }

    /**
     * public boolean isEmpty()
     * <p>
     * Returns true if and only if the heap is empty.
     */
    public boolean isEmpty() {
        return this.size()==0;
        //return this.min.getRank() == -1; //Returns true if the heap is empty
    }

    /**
     * public HeapNode insert(int key)
     * <p>
     * Creates a node (of type HeapNode) which contains the given key, and inserts it into the heap.
     * The added key is assumed not to already belong to the heap.
     * <p>
     * Returns the newly created node.
     */
    public HeapNode insert(int key) {
        HeapNode node = new HeapNode(key); //Creating a node (of type HeapNode) which contains the given key
        return insertNode(node);

    }
    public HeapNode insertNode(HeapNode node) // deal with insertion of nodes
    {

        if (this.isEmpty()|| this.min==null){// tree is empty
            this.min=node;
            this.last=node;
            this.first=node;
            node.setNext(node);
            node.setPrev(node);
            this.size++;
            this.trees++;
            return node;

        }
        if (node.getKey()<this.min.getKey())
        {
            this.min=node;
        }
        this.first.setPrev(node);
        node.setNext(this.first);
        node.setPrev(this.last);
        this.last.setNext(node);
        this.first=node;
        this.size++;
        this.trees++;
        return node;

    }
    /**
     * public void deleteMin()
     * <p>
     * Deletes the node containing the minimum key.
     */
    public void deleteMin() {
        this.size--;
        this.trees=this.trees+this.min.getRank()-1;
        if (this.min.getRank() == 0) { // If min has no child
            if (this.min==this.first && this.min==this.last)// heapNode  is empty
            {
                this.min = null;
                this.first=null;
                this.last=null;
                this.size=0;
                return;
            }
            else if (this.min==this.first)
            {
                this.first=min.getNext();
                this.first.setPrev(this.last);
                this.last.setNext(this.first);

            }
            else if (this.min==this.last)
            {
                this.last=min.getPrev();
                this.last.setNext(this.first);
                this.first.setPrev(this.last);

            }
            else {
                this.min.getPrev().setNext(min.getNext());
                this.min.getNext().setPrev(min.getPrev());
            }
            }
        else{ // We want to connect min sons to x.getPrev and x.getNext
            HeapNode child=this.min.getChild();
            HeapNode secondChild=child;
            child.setParent(null);
            do {// remove marks of child of min
                if (secondChild.getMark()==true)
                {
                    secondChild.setMark(false);
                    this.marked--;
                }
                secondChild = secondChild.getNext();

            }
            while(child!=secondChild);
            if (this.min!=this.last && this.min!=this.first) {// normal case
                child.getPrev().setNext(min.getNext());
                min.getNext().setPrev(child.getPrev());
                child.setPrev(min.getPrev());
                child.getPrev().setNext(child);

            }
            else if (this.min==this.last && this.first==this.min)// only tree is min tree
            {
                child.setParent(null);
                this.last=child.getPrev();
                this.first=child;
            }
            else if(this.min==this.first)// first tree is min
            {
                child.setParent(null);
                this.first=child;
                this.last.setNext(child);
                child.getPrev().setNext(min.getNext());
                min.getNext().setPrev(child.getPrev());
                this.first.setPrev(this.last);
            }

            else if (this.min==this.last)// last tree is min tree
            {
                child.setParent(null);
                this.last.getNext().setPrev(child.getPrev());
                this.last=child.getPrev();
                this.last.setNext(this.first);
                child.setPrev(this.min.getPrev());
                child.getPrev().setNext(child);
            }
        }

        this.succesiveLinking();

    }
    private void succesiveLinking()
    {
        int arrSize = (int)(Math.log(this.size) / Math.log(1.5)+3);;

        HeapNode [] bucket= new HeapNode[arrSize];
        for (int i = 0; i <bucket.length; i++) {
            bucket[i]=null;

        }
        this.toBuckets(bucket);
        this.fromBuckets(bucket);
    }

    private HeapNode link(HeapNode node1, HeapNode node2) {// needs to be the same rank
        this.trees--;
        totalLinks++;
        node1.setNext(node1);// initialzing pointers
        node1.setPrev(node1);
        node2.setNext(node2);
        node2.setPrev(node2);
        if (node1.getRank()==0)
        {
            if (node1.getKey()<node2.getKey())
            {
                node2.setParent(node1);
                node1.setChild(node2);
                node2.setNext(node2);
                node2.setPrev(node2);
                node1.setRank(node1.getRank()+1);
                return node1;

            }

            else
            {
                node1.setParent(node2);
                node2.setChild(node1);
                node1.setNext(node1);
                node1.setPrev(node1);
                node2.setRank(node2.getRank()+1);
                return node2;

            }
        }
         else {
            if (node1.getKey() < node2.getKey()) {
                node2.setParent(node1);
                node2.setNext(node1.getChild());
                node2.setPrev(node1.getChild().getPrev());
                node1.getChild().getPrev().setNext(node2);
                node1.getChild().setPrev(node2);
                node1.setChild(node2);
                node1.setRank(node1.getRank() + 1);
                return node1;
            }


            else
            {
                node1.setParent(node2);
                node1.setNext(node2.getChild());
                node1.setPrev(node2.getChild().getPrev());
                node2.getChild().getPrev().setNext(node1);
                node2.getChild().setPrev(node1);
                node2.setChild(node1);
                node2.setRank(node2.getRank() + 1);
                return node2;


            }

        }
        }
    private void toBuckets(HeapNode [] bucket)
    {
        HeapNode starter=this.first;
        starter.getPrev().setNext(null);
        while (starter!=null)
        {

            HeapNode y=starter;
            starter=starter.getNext();
            if (starter==this.first)
            {
                break;
            }
            while (bucket[y.getRank()]!=null)
            {
                y=link(y,bucket[y.getRank()]);

                bucket[y.getRank()-1]=null;
            }
            y.setPrev(y);
            y.setNext(y);
            y.setParent(null);
            bucket[y.getRank()]=y;
        }
    }
    private void fromBuckets(HeapNode [] bucket)
    {
        this.min=null;// initalize the heap
        this.first=null;
        this.last=null;
        this.trees=0;
        for (int i = 0; i < bucket.length; i++) {
            if (bucket[i]!=null)
            {
                this.trees++;
                if (this.min==null)
                {
                    this.min=bucket[i];
                    this.first=bucket[i];
                    this.last=bucket[i];
                    this.first.setNext(bucket[i]);
                    this.first.setPrev(bucket[i]);
                }
                else
                {
                    if (bucket[i].getKey()<this.min.getKey())
                    {
                        this.min=bucket[i];
                    }
                    bucket[i].setNext(this.first);
                    this.first.setPrev(bucket[i]);
                    bucket[i].setPrev(this.last);
                    this.last.setNext(bucket[i]);
                    this.last=bucket[i];

                }
            }

        }
    }



    /**
     * public HeapNode findMin()
     * <p>
     * Returns the node of the heap whose key is minimal, or null if the heap is empty.
     */
    public HeapNode findMin() {
        if (this.isEmpty()) { //If the heap is an empty heap
            return null;
        }
        return this.min;
    }

    /**
     * public void meld (FibonacciHeap heap2)
     * <p>
     * Melds heap2 with the current heap.
     */
    public void meld(FibonacciHeap heap2) {
        if (heap2.isEmpty())
        {
            return;
        }
        this.last.setNext(heap2.first);
        this.first.setPrev(heap2.last);
        heap2.first.setPrev(this.last);
        heap2.last.setNext(this.first);
        this.last=heap2.last;
        this.size=this.size+heap2.size;
        this.trees=this.trees+heap2.trees;
        this.marked=this.marked+heap2.marked;
        if (this.min.getKey() > heap2.min.getKey()) {
            this.min=heap2.min;
        }
    }


    /**
     * public int size()
     * <p>
     * Returns the number of elements in the heap.
     */
    public int size() {
        return this.size; // should be replaced by student code
    }

    /**
     * public int[] countersRep()
     * <p>
     * Return an array of counters. The i-th entry contains the number of trees of order i in the heap.
     * Note: The size of of the array depends on the maximum order of a tree, and an empty heap returns an empty array.
     */
    public int[] countersRep1() {
     int  count=0;
        int []array = {};
        if (this.isEmpty()){
            return array;
        }
        HeapNode start=this.first;
        Display a=new Display();
        a.display(this);
        int maxRank=0;
        do {
            if (start.getRank()>maxRank)
            {
                maxRank= start.getRank();
            }
            start=start.getNext();
            //System.out.println("stuck in a loop first one");
        }
        while(start!=this.first);

        int[] arr = new int[maxRank+1]; //initiializing the array
         start = this.first;
        do {
            arr[start.getRank()]++;
            start=start.getNext();
            System.out.println("stuck in a loop second one");

        }
        while(start!=this.first);
        return arr;
    }
    public int[] countersRep() {
        int  count=0;
        int []array = {};
        int maxRank=0;
        HeapNode start=this.first;
        if (this.isEmpty()){
            return array;
        }
        for (int i = 0; i <this.numberOfTrees();  i++) {
            if(start.getRank()>maxRank)
            {
                maxRank= start.getRank();
            }
        start=start.getNext();
        }


        int[] arr = new int[maxRank+1]; //initiializing the array
        for (int i = 0; i <this.numberOfTrees();  i++) {
            arr[start.getRank()]++;
            start = start.getNext();
        }
        return arr;
    }

    /**
     * public void delete(HeapNode x)
     * <p>
     * Deletes the node x from the heap.
     * It is assumed that x indeed belongs to the heap.
     */
    public void delete(HeapNode x) {// using deletemin and decreaseKey
        int key = x.getKey();
        int minKey = this.min.getKey();
            this.decreaseKey(x, (key - minKey) + 1);

        this.deleteMin();
    }

    /**
     * public void decreaseKey(HeapNode x, int delta)
     * <p>
     * Decreases the key of the node x by a non-negative value delta. The structure of the heap should be updated
     * to reflect this change (for example, the cascading cuts procedure should be applied if needed).
     */
    public void decreaseKey(HeapNode x, int delta) {

        if (x.getParent() == null || x.getParent().getKey() < x.getKey()-delta) {// if root does not matter
            x.setKey(x.getKey() - delta);
            if (x.getKey() < this.min.getKey()) {
                this.min = x;
            }
        }

        else {
            x.setKey(x.getKey() - delta);
            if (x.getKey() < this.min.getKey()) {
                this.min = x;
            }
            cascadingCuts(x);

        }


    }


    public void cascadingCuts(HeapNode x) {
        HeapNode parent=x.getParent();
        do {
            boolean sgn =parent.getMark();
            if (x.getMark())
            {
                this.marked--;
                x.setMark(false);
            }
            x.setParent(null);
            parent.setRank(parent.getRank() - 1);// update rank of parent
            if (x.getNext() == x) {
                parent.setChild(null);

            } else {
                parent.setChild(x.getNext());
                x.getPrev().setNext(x.getNext());
                x.getNext().setPrev(x.getPrev());
            }
            x.setNext(this.first);
            this.first.setPrev(x);
            x.setPrev(this.last);
            this.last.setNext(x);
            this.first=x;// changing x to first node
            this.trees++;
            totalCuts++;
            if (sgn==false)
            {
                if (parent.getParent()==null)
                {
                    break;
                }
                else
                {
                    parent.setMark(true);
                    this.marked++;
                    break;
                }
            }
            x=parent;
            parent=parent.getParent();

        }
        while (parent != null);

    }

    /**
     * public int potential()
     *
     * This function returns the current potential of the heap, which is:
     * Potential = #trees + 2*#marked
     *
     * In words: The potential equals to the number of trees in the heap
     * plus twice the number of marked nodes in the heap.
     */
    public int potential()
    {
        return this.trees+2*this.marked; // should be replaced by student code
    }

    /**
     * public static int totalLinks()
     *
     * This static function returns the total number of link operations made during the
     * run-time of the program. A link operation is the operation which gets as input two
     * trees of the same rank, and generates a tree of rank bigger by one, by hanging the
     * tree which has larger value in its root under the other tree.
     */
    public int numberOfTrees()
    {
        return this.trees;
    }
    public static int totalLinks()
    {
        return totalLinks; // should be replaced by student code
    }

    /**
     * public static int totalCuts()
     *
     * This static function returns the total number of cut operations made during the
     * run-time of the program. A cut operation is the operation which disconnects a subtree
     * from its parent (during decreaseKey/delete methods).
     */
    public static int totalCuts()
    {
        return totalCuts; // should be replaced by student code
    }

    /**
     * public static int[] kMin(FibonacciHeap H, int k)
     *
     * This static function returns the k smallest elements in a Fibonacci heap that contains a single tree.
     * The function should run in O(k*deg(H)). (deg(H) is the degree of the only tree in H.)
     *
     * ###CRITICAL### : you are NOT allowed to change H.
     */
     public static int[] kMin(FibonacciHeap H, int k)
    {
        int[] emptyArr = {};
        if (k == 0){
            return emptyArr;
        }
        int[] arr = new int[k];
        FibonacciHeap heap2 = new FibonacciHeap();
        arr[0] = H.min.getKey(); //Insert the min key to the array
        heap2.insert(H.min.getKey());
        heap2.findMin().setInfo(H.min);
        int counter = 1;
        HeapNode node = H.min;
        while ( counter < k){ //While the array is not full
            heap2.deleteMin();
            node = node.getChild();
            if (node != null) {
                int key = node.getKey();
                do {
                    heap2.insert(node.getKey());
                    heap2.first.setInfo(node);
                    node = node.getNext();
                }
                while (key != node.getKey());
            }
            node = heap2.findMin().getInfo();
            arr[counter] = node.getKey();
            counter ++;
        }

        return arr;
    }
    public void printLinked()
    {
     String nodes="";
        HeapNode   starter=this.first;
     while (starter!=null)
     {
         nodes=starter.toString()+"->";
         starter=starter.getNext();
     }
        System.out.println(nodes);
    }
    /**
     * public class HeapNode
     *
     * If you wish to implement classes other than FibonacciHeap
     * (for example HeapNode), do it in this file, not in another file.
     *
     */

    public static class Display { //print function
        public static final int MIN_SPACE_BETWEEN_NODES_IN_ROOT = 2;

        public static void display(FibonacciHeap heap) {
            List<StringBuilder> stringBuilders = new ArrayList<>();
            FibonacciHeap.HeapNode minNode = heap.findMin();
            buildRowDisplay(minNode, stringBuilders, 0);
            printHeap(stringBuilders);
        }

        private static void printHeap(List<StringBuilder> stringBuilders) {
            for (StringBuilder stringBuilder : stringBuilders) {
                System.out.println(stringBuilder.toString());
            }
        }

        private static void padAllToLength(List<StringBuilder> stringBuilders, int toLength) {
            for (int depth = 0; depth < stringBuilders.size(); depth++) {
                padDepth(stringBuilders, depth, toLength, true);
            }
        }

        private static void padDepth(List<StringBuilder> stringBuilders, int depth, int toLength, boolean isSpaces) {
            StringBuilder stringBuilder = getStringBuilder(stringBuilders, depth);
            int diff = toLength - stringBuilder.length();
            for (int time = diff - 1; time >= 0 ; time--) {
                char chr = isSpaces ? ' ' : time == 0 ? '>' : '-';
                stringBuilder.append(chr);
            }
        }

        private static void removeSpacesAtSuffix(List<StringBuilder> stringBuilders, int depth) {
            StringBuilder stringBuilder = getStringBuilder(stringBuilders, depth);
            int length = stringBuilder.length();
            if (stringBuilder.charAt(length -1) == ' ') {
                int startIndex = 0;
                for (int index = stringBuilder.length() - 1; index >= 0; index--) {
                    if (stringBuilder.charAt(index) != ' ') {
                        startIndex = index + 1;
                        break;
                    }
                }

                stringBuilder.replace(startIndex, length, "");
            }
        }

        private static StringBuilder getStringBuilder(List<StringBuilder> stringBuilders, int depth) {
            StringBuilder stringBuilder;
            if (depth < stringBuilders.size()) {
                stringBuilder = stringBuilders.get(depth);
            } else {
                stringBuilder = new StringBuilder();
                stringBuilders.add(depth, stringBuilder);

            }
            return stringBuilder;
        }

        private static int getMaxStringBuilderLength(List<StringBuilder> stringBuilders) {
            int maxLength = 0;
            for (StringBuilder stringBuilder :
                    stringBuilders) {
                maxLength = Math.max(maxLength, stringBuilder.length());
            }
            return maxLength;
        }

        private static void buildRowDisplay(FibonacciHeap.HeapNode node, List<StringBuilder> stringBuilders, int depth) {
            FibonacciHeap.HeapNode current = node;
            StringBuilder stringBuilder = getStringBuilder(stringBuilders, depth);
            do {
                buildNodeDisplay(current, stringBuilders, depth);

                removeSpacesAtSuffix(stringBuilders, depth);
                int padToLength = Math.max(stringBuilder.length() + MIN_SPACE_BETWEEN_NODES_IN_ROOT, getMaxStringBuilderLength(stringBuilders));
                padDepth(stringBuilders, depth, padToLength, current.next == node);

                padAllToLength(stringBuilders, stringBuilder.length());
                current = current.next;
            } while (current != node);
        }


        private static void buildNodeDisplay(FibonacciHeap.HeapNode node, List<StringBuilder> stringBuilders, int depth) {
            StringBuilder stringBuilder = getStringBuilder(stringBuilders, depth);
            stringBuilder.append(node.key);
            if (node.child != null) {
                StringBuilder lineSep = getStringBuilder(stringBuilders, depth + 1);
                padDepth(stringBuilders, depth + 1, stringBuilder.length() - 1, true);
                padDepth(stringBuilders, depth + 2, stringBuilder.length() - 1, true);
                lineSep.append('|');
                buildRowDisplay(node.child, stringBuilders, depth + 2);
            }
        }
    }
    public boolean checkCounters()
    {
        int [] arr=this.countersRep();
        HeapNode start=this.first;
        do {
            int rank=start.getRank();
            int realRank=0;
            HeapNode child=start.getChild();
            if (child==null)
            {
                arr[realRank]--;
                start=start.getNext();
                continue;
            }
            child.getPrev().setNext(null);
            while(child!=null)
            {
                realRank++;
                child=child.getNext();
            }
            start.getChild().getPrev().setNext(start.getChild());
            arr[realRank]--;
            start=start.getNext();
        }
        while(start!=this.first);
        for (int i = 0; i <arr.length ; i++) {
            if (arr[i]!=0)
            {
                return false;
            }

        }
        return true;
    }
    public boolean checkRank()// checks rank of tree is correct
    {

        HeapNode start=this.first;
        do {
           int rank=start.getRank();
           int realRank=0;
           HeapNode child=start.getChild();
          if (child==null)
          {
              if (rank!=0)
              {
                  System.out.println("rank is not correct");
                  return false;
              }
              continue;
          }
          child.getPrev().setNext(null);
          while(child!=null)
           {
               realRank++;
               child=child.getNext();
           }
          start.getChild().getPrev().setNext(start.getChild());
            if(realRank!=rank)
            {
                System.out.println("rank is not correct");
                return false;
            }
            start=start.getNext();
        }
        while(start!=this.first);

        return true;
    }
    public boolean checkParent()
    {
        HeapNode start = this.first;
        do {


            if (start.getParent() != null) {
                System.out.println("parent is not null");
                return false;
            }
        start=start.getNext();
        }
        while(start!=this.first);
        return true;
    }
    public static class HeapNode{

        private int key;
        private int rank;
        private boolean mark;
        private HeapNode child;
        private HeapNode next;
        private HeapNode prev;
        private HeapNode parent;
        private HeapNode info;
        private HeapNode() {
            this.key = Integer.MAX_VALUE;
            this.rank = -1;
            this.mark = false;
            this.child = null;
            this.next = null;
            this.prev = null;
            this.parent = null;
            this.info=null;
        }
        public HeapNode(int key) {
            this.key = key;
            this.rank = 0;
            this.mark = false;
            this.child = null;
            this.next = null;
            this.prev = null;
            this.parent = null;
            this.info=null;

        }

        public int getKey() {
            return this.key;
        }
        public int getRank(){
            return this.rank;
        }
        public boolean getMark(){
            return this.mark;
        }
        public HeapNode getChild(){
            return this.child;
        }
        public HeapNode getNext(){
            return this.next;
        }
        public HeapNode getPrev(){
            return this.prev;
        }
        public HeapNode getParent(){
            return this.parent;
        }
        public HeapNode getInfo() {return this.info;}
        public void setKey(int key){
            this.key = key;
        }
        public void setRank(int rank){
            this.rank = rank;
        }
        public void setMark(boolean mark){
            this.mark = mark;
        }
        public void setChild(HeapNode child){
            this.child = child;
        }
        public void setNext(HeapNode next){
            this.next = next;
        }
        public void setPrev(HeapNode prev){
            this.prev = prev;
        }
        public void setParent(HeapNode parent){
            this.parent = parent;
        }
        public void setInfo(HeapNode info)
        {
            this.info=info;
    }
        public String toString(){
            return ""+this.key;
        }
    }
}