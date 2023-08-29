template< typename T > struct myallocator : public std::allocator<T>
{
  typedef std::allocator<T> base ;
  typedef typename base::size_type size_type;
  typedef typename base::difference_type  difference_type;
  typedef typename base::pointer pointer;
  typedef typename base::const_pointer const_pointer;
  typedef typename base::reference reference;
  typedef typename base::const_reference const_reference;
  typedef typename base::value_type value_type;
  myallocator() throw() {}
  myallocator( const myallocator& a ) throw() : base(a) {}
  template <typename X> myallocator(
             const myallocator<X>& ) throw() {}
  ~myallocator() throw() {}

  template <typename X> struct rebind
  { typedef myallocator<X> other; };

  pointer allocate( size_type sz, const void* hint = 0 )
  {
    // record alloc request eg. ++num_allocs ;
    return base::allocate(sz,hint) ;
  }
  void deallocate( pointer p, size_type n )
  {
    // record dealloc request eg. --num_allocs ;
    return base::deallocate(p,n) ;
  }
};

template<typename T> inline bool operator==(
               const myallocator<T>&, const myallocator<T>& )
{ return true; }

template<typename T> inline bool operator!=(
               const myallocator<T>&, const myallocator<T>& )
{ return false; }

