//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of pSpatiocyte
//
//        Copyright (C) 2019 Satya N.V. Arjunan
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Motocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Motocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Motocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// based on DynamicPriorityQueue.hpp of E-Cell Project
// by Koichi Takahashi, Eiichiro Adachi and Moriyoshi Koizumi
//


#ifndef __PriorityQueue_hpp
#define __PriorityQueue_hpp

class VolatileIDPolicy {
public:
  typedef size_t    Index;
  typedef Index     ID;

  class IDIterator {
  public:
    typedef ptrdiff_t difference_type;
    typedef VolatileIDPolicy::Index size_type;
    typedef VolatileIDPolicy::ID value_type; 
  public:
    IDIterator( VolatileIDPolicy::Index _idx ):
      idx( _idx ) {} 
    IDIterator operator++(int) {
      IDIterator retval = *this;
      ++this->idx;
      return retval;
    }
    IDIterator& operator++() {
      ++this->idx;
      return *this;
    } 
    IDIterator operator--(int) {
      IDIterator retval = *this;
      --this->idx;
      return retval;
    }
    IDIterator& operator--() {
      --this->idx;
      return *this;
    }
    value_type operator*() {
      return this->idx;
    }
    difference_type operator-( const IDIterator& that ) {
      return this->idx - that.idx;
    }
    bool operator==( const IDIterator& that ) {
      return this->idx == that.idx;
    }
    bool operator!=( const IDIterator& that ) {
      return this->idx != that.idx;
    }
    bool operator <( const IDIterator& that ) {
      return this->idx < that.idx;
    }
    bool operator >( const IDIterator& that ) {
      return this->idx > that.idx;
    }
  private:
    size_type idx;
  };
public:
  void reset() {}
  void clear() {}
  Index getIndex( ID id ) const {
    return id;
  }
  ID getIDByIndex( Index index ) const {
    return index;
  }
  ID push( Index index ) {
    maxIndex = index + 1;
    return index;
  }
  void pop( Index index ) {
    maxIndex = index;
  }
  bool checkConsistency( Index size ) const {
    return true;
  }
  IDIterator begin() const {
    return IDIterator( 0 );
  }
  IDIterator end() const {
    return IDIterator( maxIndex );
  }
private:
  Index maxIndex;
};

template<typename Item, class IDPolicy = VolatileIDPolicy>
class PriorityQueue {
public:
  typedef std::vector<Item> ItemVector;
  typedef typename IDPolicy::ID ID;
  typedef typename IDPolicy::Index Index;
  typedef std::vector<Index> IndexVector;
  typedef typename IDPolicy::IDIterator IDIterator;
public:
  void clear(); 
  inline ID push(const Item& item); 
  inline void movePos(Index pos); 
  PriorityQueue() {}
  bool isEmpty() const {
    return this->itemVector.empty();
  } 
  Index getSize() const {
    return this->itemVector.size();
  }
  Item& get_top() {
    return this->itemVector[getTopIndex()];
  } 
  Item const& get_top() const {
    return this->itemVector[getTopIndex()];
  } 
  Item& get(const ID id) {
    return this->itemVector[this->pol.getIndex(id)];
  } 
  Item const& get(const ID id) const {
    return this->itemVector[this->pol.getIndex(id)];
  } 
  ID get_top_id() const {
    return this->pol.getIDByIndex(getTopIndex());
  } 
  Item& operator[](const ID id) {
    return get(id);
  } 
  Item const& operator[](const ID id) const {
    return get(id);
  }
  Item& getByIndex(const Index index) {
    return this->itemVector[index];
  } 
  Index getTopIndex() const {
    return this->heap[0];
  } 
  void move(Index index) {
    const Index pos(this->positionVector[index]);
    movePos(pos);
  } 
  void move_top() {
    move_down_pos(0);
  } 
  void moveUpByIndex(Index index) {
    const Index position(this->positionVector[index]);
    move_up_pos(position);
  } 
  void move_up(ID id) {
    moveUpByIndex(pol.getIndex(id));
  } 
  void moveDownByIndex(Index index) {
    const Index position(this->positionVector[index]);
    move_down_pos(position);
  } 
  void move_down(ID id) {
    moveDownByIndex(pol.getIndex(id));
  } 
  IDIterator begin() const {
    return pol.begin();
  } 
  IDIterator end() const {
    return pol.end();
  }
private:
  inline void move_up_pos(Index position, Index start = 0);
  inline void move_down_pos(Index position); 
private:
  ItemVector itemVector;
  IndexVector heap; 
  IndexVector positionVector;
  IDPolicy pol;
};

template<typename Item, class IDPolicy>
void PriorityQueue<Item, IDPolicy>::clear() { 
  this->itemVector.clear();
  this->heap.clear();
  this->positionVector.clear();
  this->pol.clear();
}

template<typename Item, class IDPolicy>
void PriorityQueue<Item, IDPolicy>::movePos(Index pos) {
  const Index index(this->heap[pos]);
  const Item& item(this->itemVector[index]); 
  const Index size(getSize()); 
  if(pos < size/2) {
    const Index succ(2*pos+1);
    if(this->itemVector[this->heap[succ]]->get_time() < item->get_time() ||
       (this->itemVector[this->heap[succ]]->get_time() == item->get_time() &&
        this->itemVector[this->heap[succ]]->get_priority() >
        item->get_priority()) ||
       (succ+1 < size &&
        this->itemVector[this->heap[succ+1]]->get_time() < item->get_time()) ||
       (succ+1 < size &&
        this->itemVector[this->heap[succ+1]]->get_time() == item->get_time() &&
        this->itemVector[this->heap[succ+1]]->get_priority() >
        item->get_priority())) {
      move_down_pos(pos);
      return;
    }
  } 
  if(pos > 0) {
    const Index pred((pos-1)/2);
    if(item->get_time() < this->itemVector[this->heap[pred]]->get_time() ||
       (item->get_time() == this->itemVector[this->heap[pred]]->get_time() &&
        item->get_priority() >
        this->itemVector[this->heap[pred]]->get_priority())) {
      move_up_pos(pos);
      return;
    }
  }
}

template<typename Item, class IDPolicy>
void PriorityQueue<Item, IDPolicy>::move_up_pos(Index position, Index start) {
  if(position == 0) {
    return;
  } 
  const Index index(this->heap[position]);
  const Item& item(this->itemVector[index]); 
  Index pos(position);
  while(pos > start) {
    const Index pred((pos-1)/2);
    const Index predIndex(this->heap[pred]);
    if(this->itemVector[predIndex]->get_time() < item->get_time() ||
       (this->itemVector[predIndex]->get_time() == item->get_time() &&
        this->itemVector[predIndex]->get_priority() > item->get_priority())) {
      break;
    } 
    this->heap[pos] = predIndex;
    this->positionVector[predIndex] = pos;
    pos = pred;
  } 
  this->heap[pos] = index;
  this->positionVector[index] = pos;
}

template<typename Item, class IDPolicy>
void PriorityQueue<Item, IDPolicy>::move_down_pos(Index position) {
  const Index index(this->heap[position]);
  const Index size(getSize()); 
  Index succ(2*position+1);
  Index pos(position);
  while(succ < size) {
    const Index rightPos(succ + 1);
    if(rightPos < size &&
       (this->itemVector[this->heap[rightPos]]->get_time() <
        this->itemVector[this->heap[succ]]->get_time() ||
        (this->itemVector[this->heap[rightPos]]->get_time() ==
         this->itemVector[this->heap[succ]]->get_time() &&
         this->itemVector[this->heap[rightPos]]->get_priority() >
         this->itemVector[this->heap[succ]]->get_priority()))) {
      succ = rightPos;
    } 
    this->heap[pos] = this->heap[succ];
    this->positionVector[this->heap[pos]] = pos;
    pos = succ;
    succ = 2*pos+1;
  } 
  this->heap[pos] = index;
  this->positionVector[index] = pos; 
  move_up_pos(pos, position);
} 

template<typename Item, class IDPolicy> 
typename PriorityQueue<Item, IDPolicy>::ID
PriorityQueue<Item, IDPolicy>::push(const Item& item) {
  const Index index(getSize()); 
  this->itemVector.push_back(item);
  this->heap.push_back(index);
  this->positionVector.push_back(index); 
  const ID id(this->pol.push(index));
  move_up_pos(index); 
  return id;
}

#endif // __PriorityQueue_hpp
