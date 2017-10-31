/*
 * CuckooHashSetTest.java
 * JUnit based test
 *
 * Created on May 12, 2006, 11:30 PM
 */

package hash;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import junit.framework.*;
import java.util.List;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author mph
 */
public class RefinableCuckooHashSetTest extends TestCase {
  private final static int THREADS = 8;
  private final static int TEST_SIZE = 64;
  private final static int PER_THREAD = TEST_SIZE / THREADS;
  int index;
  RefinableCuckooHashSet<Integer> instance;
  Thread[] thread = new Thread[THREADS];
  Random random = new Random(0);
  AtomicInteger missed = new AtomicInteger(0);
  
  public RefinableCuckooHashSetTest(String testName) {
    super(testName);
    instance = new RefinableCuckooHashSet<Integer>(8);
  }
  
  public static Test suite() {
    TestSuite suite = new TestSuite(RefinableCuckooHashSetTest.class);
    return suite;
  }
  
  /**
   * Sequential calls.
   */
  public void testSequential() {
    System.out.println("sequential add, contains, and remove");
    Set<Integer> log = new HashSet<Integer>(TEST_SIZE);

    for (int i = 0; i < TEST_SIZE; i++) {
      int x = random.nextInt();
      if (instance.add(x)) {
        log.add(x);
      }
    }
    if (!instance.check()) {
      fail("sequential: instance inconsistent");
    }
    for (int i : log) {
      if (!instance.contains(i)) {
        fail("sequential: expected to find: " + i );
      }
    }

    for (int i : log) {
      if (!instance.remove(i)) {
        fail("sequential: bad remove: " + i );
      }
    }
  }
  /**
   * Parallel adds, sequential removes
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
    if (!instance.check(TEST_SIZE)) {
      fail("concurrent add: instance inconsistent");
    }
    for (int i = 0; i < TEST_SIZE; i++) {
      if (!instance.contains(i)) {
        fail("contains expected to find: " + i );
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
    if (!instance.check(0)) {
      fail("concurrent remove: instance inconsistent");
    }
  }
  /**
   * Sequential enqueues, parallel dequeues
   */
  public void testParallelBoth()  throws Exception {
    System.out.println("parallel both");
    Thread[] myThreads = new Thread[2 * THREADS];
    for (int i = 0; i < THREADS; i++) {
      myThreads[i] = new AddThread(i * PER_THREAD);
      myThreads[i + THREADS] = new RemoveThread2(i * PER_THREAD);
    }
    for (int i = 0; i < 2 * THREADS; i ++) {
      myThreads[i].start();
    }
    for (int i = 0; i < 2 * THREADS; i ++) {
      myThreads[i].join();
    }
    if (!instance.check(missed.get())) {
      fail("concurrent both: instance inconsistent");
    }
  }
  class AddThread extends Thread {
    int value;
    AddThread(int i) {
      value = i;
    }
    public void run() {
      for (int i = 0; i < PER_THREAD; i++) {
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
        if (!instance.remove(value+i)) {
          fail("Th " + value + "\tDeqThread: missing value: " + (value + i));
        }
      }
    }
  }
  class RemoveThread2 extends Thread {
    int value;
    RemoveThread2(int i) {
      value = i;
    }
    public void run() {
      for (int i = 0; i < PER_THREAD; i++) {
        if (!instance.remove(value+i)) {
          missed.getAndIncrement();
        }
      }
    }
  }
}
