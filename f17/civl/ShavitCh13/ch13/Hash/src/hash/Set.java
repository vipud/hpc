/*
 * Set.java
 *
 * Created on December 29, 2005, 11:43 PM
 *
 * From "The Art of Multiprocessor Programming",
 * by Maurice Herlihy and Nir Shavit.
 * Copyright 2006 Elsevier Inc. All rights reserved.
 */

package hash;

import java.util.Iterator;
/**
 * Interface satisfied by various buckets
 * @author Maurice Herlihy
 */
public interface Set<T> extends Iterable<T> {
  /**
   * add object with given key
   * @param x object to add
   * @return whether object was absent
   */
  boolean add(T x);
  /**
   * remove object from bucket
   * @param x object to remove
   * @return whether object was found
   */
  public boolean remove(T x);
  /**
   * is object in bucket?
   * @param x object being sought
   * @return whether object is present
   */
  public boolean contains(T x);
  /**
   * iterate through objects in bucket
   * @return iterator over elements
   */
  public Iterator<T> iterator();
}
