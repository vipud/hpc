/*
 * StripedHashSet.java
 *
 * Created on November 15, 2006, 3:35 PM
 *
 * From "The Art of Multiprocessor Programming",
 * by Maurice Herlihy and Nir Shavit.
 * Copyright 2006 Elsevier Inc. All rights reserved.
 */

package hash;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Concurrent Cuckoo hashing using lock striping.
 * @param T Item type.
 * @author Maurice Herlihy
 */
public class StripedHashSet<T> extends BaseHashSet<T>{
  final Lock[] locks;
  
  public StripedHashSet(int capacity) {
    super(capacity);
    locks  = new Lock[capacity];
    for (int j = 0; j < locks.length; j++) {
      locks[j] = new ReentrantLock();
    }
  }
  
  /**
   * double the set size
   */
  public void resize() {
    int oldCapacity = table.length;
    for (Lock lock : locks) {
      lock.lock();
    }
    try {
      if (oldCapacity != table.length) {
        return; // someone beat us to it
      }
      int newCapacity  = 2 * oldCapacity;
      List<T>[] oldTable = table;
      table = (List<T>[]) new List[newCapacity];
      for (int i = 0; i < newCapacity; i++)
        table[i] = new ArrayList<T>();
      initializeFrom(oldTable);
    } finally {
      for (Lock lock : locks) {
        lock.unlock();
      }
    }
  }
  private void initializeFrom(List<T>[] oldTable) {
    for (List<T> bucket : oldTable) {
      for (T x : bucket) {
        int myBucket = Math.abs(x.hashCode() % table.length);
        table[myBucket].add(x);
      }
    }
  }
  /**
   * Synchronize before adding, removing, or testing for item
   * @param x item involved
   */
  public final void acquire(T x) {
    int myBucket = Math.abs(x.hashCode() % locks.length);
    locks[myBucket].lock();
  }
  /**
   * synchronize after adding, removing, or testing for item
   * @param x item involved
   */
  public void release(T x) {
    int myBucket = Math.abs(x.hashCode() % locks.length);
    locks[myBucket].unlock();
  }
  public boolean policy() {
    return size / table.length > 4;
  }
}
