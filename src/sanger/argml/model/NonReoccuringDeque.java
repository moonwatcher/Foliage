/* 
 * Foliage. An Ancestral Recombination Graph Manipulation Library.
 * 
 * Copyright (c) 2008 Genome Research Ltd.
 * 
 * Author: Lior Galanti <lior.galanti@gmail.com>
 * 
 * This file is part of Foliage.
 * Foliage is free software; you can redistribute it and/or modify it under the terms of 
 * the GNU General Public License as published by the Free Software Foundation; 
 * either version 2 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 * 
 */

package sanger.argml.model;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;

public class NonReoccuringDeque<E> implements Deque<E> {

	private HashSet<E> occured;
	private ArrayDeque<E> deque;
	
	public NonReoccuringDeque(){
		this.occured = new HashSet<E>();
		this.deque = new ArrayDeque<E>();
	}
	
	public boolean add(E e) {
		return this.occured.add(e) ? this.deque.add(e) : false;
	}

	public void addFirst(E e) {
		if ( this.occured.add(e) ) { this.deque.addFirst(e); }
	}

	public void addLast(E e) {
		if ( this.occured.add(e) ) { this.deque.addLast(e); }
	}

	public boolean contains(Object o) {
		return this.occured.contains(o);
	}

	public Iterator<E> descendingIterator() {
		return this.deque.descendingIterator();
	}

	public E element() {
		return this.deque.element();
	}

	public E getFirst() {
		return this.deque.getFirst();
	}

	public E getLast() {
		return this.deque.getLast();
	}

	public Iterator<E> iterator() {
		return this.deque.iterator();
	}

	public boolean offer(E e) {
		boolean result = this.occured.add(e); 
		if(result){
			result = this.deque.offer(e);
			if(!result){this.occured.remove(e);}
		}
		return result;
	}

	public boolean offerFirst(E e) {
		boolean result = this.occured.add(e); 
		if(result){
			result = this.deque.offerFirst(e);
			if(!result){this.occured.remove(e);}
		}
		return result;
	}

	public boolean offerLast(E e) {
		boolean result = this.occured.add(e); 
		if(result){
			result = this.deque.offerLast(e);
			if(!result){this.occured.remove(e);}
		}
		return result;
	}

	public E peek() {
		return this.deque.peek();
	}

	public E peekFirst() {
		return this.deque.peekFirst();
	}

	public E peekLast() {
		return this.deque.peekLast();
	}

	public E poll() {
		return this.deque.poll();
	}

	public E pollFirst() {
		return this.deque.pollFirst();
	}

	public E pollLast() {
		return this.deque.pollLast();
	}

	public E pop() {
		return this.deque.pop();
	}

	public void push(E e) {
		if ( this.occured.add(e) ) { this.deque.push(e);} 
	}

	public E remove() {
		return this.deque.remove();
	}

	public boolean remove(Object o) {
		return this.deque.remove(o);
	}

	public E removeFirst() {
		return this.deque.removeFirst();
	}

	public boolean removeFirstOccurrence(Object o) {
		return this.deque.removeFirstOccurrence(o);
	}

	public E removeLast() {
		return this.deque.removeLast();
	}

	public boolean removeLastOccurrence(Object o) {
		return this.deque.removeLastOccurrence(o);
	}

	public int size() {
		return this.deque.size();
	}

	public boolean addAll(Collection c) {
		return false;
	}

	public void clear() {
		this.deque.clear();
		this.occured.clear();
	}

	public boolean containsAll(Collection c) {
		return this.deque.containsAll(c);
	}

	public boolean isEmpty() {
		return this.deque.isEmpty();
	}

	public boolean removeAll(Collection c) {
		return this.deque.retainAll(c);
	}

	public boolean retainAll(Collection c) {
		return this.deque.retainAll(c);
	}

	public Object[] toArray() {
		return this.deque.toArray();
	}

	public <T> T[] toArray(T[] a) {
		return this.deque.toArray(a);
	}

	public String toString(){
		return deque.toString();
	}
}
