/*
 * CoarseHashSet.java
 *
 * Created on December 29, 2005, 11:50 PM
 *
 * From "The Art of Multiprocessor Programming",
 * by Maurice Herlihy and Nir Shavit.
 * Copyright 2006 Elsevier Inc. All rights reserved.
 */
package hash;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
public class ListBucket<T> implements Set<T> {
  private List<T> list = new LinkedList<T>();
  // add object with given key
  public boolean add(T x) {
    return list.add(x);
  }
  // remove object with given key
  public boolean remove(T x) {
    return list.remove(x);
  }
  // is object present?
  public boolean contains(T x) {
    return list.contains(x);
  }
  // iterate over Set elements
  public Iterator<T> iterator() {
    return list.iterator();
  }
}
