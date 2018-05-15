/*
 * FineHashSetTest.java
 * JUnit based test
 *
 * Created on December 30, 2005, 12:14 AM
 */
package lists;

import junit.framework.*;

/**
 * @author Maurice Herlihy
 */
public class LockFreeListTest extends TestCase {
  
  private final static int THREADS = 8;
  private final static int TEST_SIZE = 128;
  private final static int PER_THREAD = TEST_SIZE / THREADS;
  LockFreeList<Integer> instance;
  Thread[] thread = new Thread[THREADS];
  
  public LockFreeListTest(String testName) {
    super(testName);
    instance = new LockFreeList<Integer>();
  }
  
  public static Test suite() {
    TestSuite suite = new TestSuite(LockFreeListTest.class);
    
    return suite;
  }
  
  /**
   * Sequential calls.
   */
  public void testSequential() {
    System.out.println("sequential add, contains, and remove");
    
    for (int i = 0; i < TEST_SIZE; i++) {
      instance.add(i);
    }
    for (int i = 0; i < TEST_SIZE; i++) {
      if (!instance.contains(i)) {
        fail("bad contains: " + i );
      }
    }
    for (int i = 0; i < TEST_SIZE; i++) {
      if (!instance.remove(i)) {
        fail("bad remove: " + i );
      }
    }
  }
  /**
   * Parallel add, sequential removes
   * @throws java.lang.Exception 
   */
  public void testParallelAdd()  throws Exception {
    System.out.println("parallel add");
    for (int i = 0; i < THREADS; i++) {
      thread[i] = new AddThread(i * PER_THREAD);
    }
    for (int i = 0; i < THREADS; i ++) {
      thread[i].start();
    }
    for (int i = 0; i < THREADS; i ++) {
      thread[i].join();
    }
    for (int i = 0; i < TEST_SIZE; i++) {
      if (!instance.contains(i)) {
        fail("bad contains: " + i );
      }
    }
    for (int i = 0; i < TEST_SIZE; i++) {
      if (!instance.remove(i)) {
        fail("bad remove: " + i );
      }
    }
  }

  /**
   * Sequential adds, parallel removes
   * @throws java.lang.Exception 
   */
  public void testParallelRemove()  throws Exception {
    System.out.println("parallel remove");
    for (int i = 0; i < TEST_SIZE; i++) {
      instance.add(i);
    }
    for (int i = 0; i < TEST_SIZE; i++) {
      if (!instance.contains(i)) {
        fail("bad contains: " + i );
      }
    }
    for (int i = 0; i < THREADS; i++) {
      thread[i] = new RemoveThread(i * PER_THREAD);
    }
    for (int i = 0; i < THREADS; i ++) {
      thread[i].start();
    }
    for (int i = 0; i < THREADS; i ++) {
      thread[i].join();
    }
  }
  
  /**
   * Parallel adds, removes
   * @throws java.lang.Exception 
   */
  public void testParallelBoth()  throws Exception {
    System.out.println("parallel both");
    Thread[] myThreads = new Thread[2 * THREADS];
    for (int i = 0; i < THREADS; i++) {
      myThreads[i] = new AddThread(i * PER_THREAD);
      myThreads[i + THREADS] = new RemoveThread(i * PER_THREAD);
    }
    for (int i = 0; i < 2 * THREADS; i ++) {
      myThreads[i].start();
    }
    for (int i = 0; i < 2 * THREADS; i ++) {
      myThreads[i].join();
    }
  }
  class AddThread extends Thread {
    int value;
    AddThread(int i) {
      value = i;
    }
    public void run() {
      for (int i = 0; i < PER_THREAD; i++) {
        if (value + i == 128) {
          int z =0;
        }
        instance.add(value + i);
      }
    }
  }
  class RemoveThread extends Thread {
    int value;
    RemoveThread(int i) {
      value = i;
    }
    public void run() {
      for (int i = 0; i < PER_THREAD; i++) {
        if (!instance.remove(value + i)) {
          fail("RemoveThread: duplicate remove: " + (value + i));
        }
      }
    }
  }
  
}
