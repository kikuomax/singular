#ifndef _SINGULAR_VECTOR_H
#define _SINGULAR_VECTOR_H

#include "singular/Exception.h"

#include <cstddef>
#include <iterator>
#include <sstream>

namespace singular {

	/** Vector. */
	template < typename T >
	class Vector {
	private:
		/** Memory block of this vector. */
		T* pBlock;

		/** Size of this vector. */
		size_t len;

		/** Distance to jump from one element to the next element. */
		ptrdiff_t delta;
	public:
		/** Const iterator. */
		class const_iterator
			: public std::iterator< std::bidirectional_iterator_tag, const T >
		{
		private:
			friend class Vector;

			/** Pointer to an element. */
			const T* ptr;

			/** Distance to jump to the next element. */
			ptrdiff_t delta;
		private:
			/**
			 * Constructs an iterator with given pointer and delta.
			 *
			 * @param ptr
			 *     Pointer to an element.
			 * @param delta
			 *     Distance to jump to the next element.
			 */
			inline const_iterator(const T* ptr, ptrdiff_t delta)
				: ptr(ptr), delta(delta) {}

			/**
			 * Returns an iterator that points to the first element in a given
			 * vector.
			 *
			 * @param v
			 *     Vector in which elements are to be iterated.
			 * @return
			 *     Iterator that points to the first element in `v`.
			 */
			static inline const_iterator begin(const Vector< T >& v) {
				return const_iterator(v.pBlock, v.delta);
			}

			/**
			 * Returns an iterator that points the stop element in a given
			 * vector.
			 *
			 * @param v
			 *     Vector in which elements are to be iterated.
			 * @return
			 *     Stop iterator of `pVector`.
			 */
			static inline const_iterator end(const Vector< T >& v) {
				const T* ptr = v.pBlock + v.size() * v.delta;
				return const_iterator(ptr, v.delta);
			}

			/**
			 * Moves the pointer.
			 *
			 * @param delta
			 *     Amount of the move.
			 */
			inline void move(ptrdiff_t delta) {
				this->ptr += delta;
			}
		public:
			/**
			 * Returns the element that this iterator points to.
			 *
			 * @return
			 *     Element that this iterator points to.
			 */
			inline const T& operator *() const {
				return *this->ptr;
			}

			/**
			 * Returns the pointer of this iterator.
			 *
			 * @return
			 *     Pointer of this iterator.
			 */
			inline const T* operator ->() const {
				return this->ptr;
			}

			/**
			 * Moves forward this iterator (pre-increment).
			 *
			 * @return
			 *     This iterator that points to the next element.
			 */
			inline const_iterator& operator ++() {
				this->move(this->delta);
				return *this;
			}

			/**
			 * Moves forward this iterator (post-increment).
			 *
			 * @return
			 *     Iterator that points to the element before the move.
			 */
			const_iterator operator ++(int) {
				const_iterator prev = *this;
				this->move(this->delta);
				return prev;
			}

			/**
			 * Moves backward this iterator (pre-decrement).
			 *
			 * @return
			 *     This iterator that points to the previous element.
			 */
			inline const_iterator& operator --() {
				this->move(-this->delta);
				return *this;
			}

			/**
			 * Moves backward this iterator (post-decrement).
			 *
			 * @return
			 *     Iterator that points to the element before the move.
			 */
			const_iterator operator --(int) {
				const_iterator next = *this;
				this->move(-this->delta);
				return next;
			}

			/**
			 * Returns whether this iterator and a given iterator point
			 * the same element.
			 *
			 * @param rhs
			 *     Iterator to be tested.
			 * @return
			 *     Whether this iterator and `rhs` point the same element.
			 */
			inline bool operator ==(const const_iterator& rhs) const {
				return this->ptr == rhs.ptr;
			}

			/**
			 * Returns whether this iterator and a given iterator point
			 * different elements.
			 *
			 * @param rhs
			 *     Iterator to be tested.
			 * @return
			 *     Whether this iterator and `rhs` point different elements.
			 */
			inline bool operator !=(const const_iterator& rhs) const {
				return this->ptr != rhs.ptr;
			}
		};
	public:
		/**
		 * Constructs an instance that wraps a given memory block.
		 *
		 * A vector never holds the ownership of `pBlock`.
		 *
		 * @param pBlock
		 *     Pointer to a memory block of the vector.
		 * @param size
		 *     Size of the vector.
		 * @param delta
		 *     Distance to jump from one element to the next element.
		 */
		Vector(T* pBlock, size_t size, ptrdiff_t delta)
			: pBlock(pBlock), len(size), delta(delta) {}

		/**
		 * Returns the size of this vector.
		 *
		 * @return
		 *     Size of this vector.
		 */
		inline size_t size() const {
			return this->len;
		}

		/**
		 * Returns the element at a given index in this vector.
		 *
		 * Undefined if `idx >= this.size()`.
		 *
		 * @param idx
		 *     Index of the element to be obtained.
		 * @return
		 *     Element at `idx`.
		 *     Changes to a returned element is reflected to this vector.
		 */
		inline T& operator [](size_t idx) {
			return this->pBlock[idx * this->delta];
		}

		/**
		 * Returns the element at a given index in this vector.
		 *
		 * Undefined if `idx >= this.size()`.
		 *
		 * @param idx
		 *     Index of the element to be obtained.
		 * @return
		 *     Element at `idx`.
		 */
		inline const T& operator [](size_t idx) const {
			return this->pBlock[idx * this->delta];
		}

		/**
		 * Starts iterating on elements in this vector.
		 *
		 * @return
		 *     Iterator that iterates elements in this vector.
		 */
		const_iterator begin() const {
			return const_iterator::begin(*this);
		}

		/**
		 * Returns the stop iterator of this vector.
		 *
		 * @return
		 *     Stop iterator of this vector.
		 */
		const_iterator end() const {
			return const_iterator::end(*this);
		}

		/**
		 * Returns a subvector of this vector.
		 *
		 * @param start
		 *     Index of the beginning of the slice.
		 * @return
		 *     Sliced vector.
		 * @throws Exception
		 *     If `start > this->size()`.
		 */
		Vector< T > slice(size_t start) const {
			if (this->size() < start) {
				std::ostringstream msg;
				msg << "start must be <= this->size() but "
					<< start << " > " << this->size();
				throw Exception(msg.str());
			}
			return Vector< T >(this->pBlock + start * this->delta,
							   this->size() - start,
							   this->delta);
		}

		/**
		 * Implicit conversion to a const variant.
		 *
		 * @return
		 *     Const variant of this vector.
		 */
		operator Vector< const T >() const {
			return Vector< const T >(this->pBlock, this->size(), this->delta);
		}
	};

}

#endif
