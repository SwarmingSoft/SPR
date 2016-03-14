#ifndef LIBMEMORY_HEADER
#define LIBMEMORY_HEADER

#include <algorithm> //copy
#include <cassert> //assert
#include <cstring> //memset

template <typename X, typename Y> auto intmod(X x, Y y) -> decltype(((x %= y) < 0) ? x+y : x)
{
    assert(y >= 0);
    return ((x %= y) < 0) ? x+y : x;
    //return ((x % y) + y) % y; //slower by a factor of 2
}

namespace MemoryConst
{
    static const int Dynamic = -1;
}

template <class T> struct TCOORDINATE
{
    T X, Y;

    TCOORDINATE(void): X(), Y() { }

    TCOORDINATE(T X, T Y) : X(X), Y(Y) { }

    TCOORDINATE(const TCOORDINATE<T> &rhs) : X(rhs.X), Y(rhs.Y) {}

    bool operator== (const TCOORDINATE<T> &rhs)
    {
        return (X == rhs.X && Y == rhs.Y);
    }

    virtual ~TCOORDINATE(void) { }
};

typedef TCOORDINATE<unsigned int> COORDINATE;


template <class T, int BufferSize> class Memory
{
protected:

    T Buffer[BufferSize];

public:

    Memory(void) {}
    Memory(unsigned int) {} //?? Ignore Value

    virtual ~Memory(void) {}

    bool operator==(const Memory &rhs) const
    {
        return std::equal(Buffer, Buffer + BufferSize, rhs.Buffer);
    }

    inline unsigned int Size(void)
    {
        return BufferSize;
    }

    inline const T *data(void) const
    {
        return Buffer;
    }

    inline T *data(void)
    {
        return Buffer;
    }

    //DEBUG
    inline bool IsDynamic(void) const
    {
        return false;
    }
};

template <class T> class Memory<T, 0>
{
public:
    Memory(void) {}
    Memory(unsigned int) {} //?? Ignore Value

    virtual ~Memory(void) {}

    bool operator==(const Memory &rhs) const
    {
        return true;
    }

    inline unsigned int Size(void)
    {
        return 0;
    }

    inline const T *data(void) const
    {
        return 0;
    }

    inline T *data(void)
    {
        return 0;
    }

    //DEBUG
    inline bool IsDynamic(void) const
    {
        return false;
    }
    //DEBUG
};

template <class T> class Memory<T, MemoryConst::Dynamic>
{
protected:

    T *Buffer;
    unsigned int size;

public:
//private:

    //1: does not work somehow, Memory() can still be created
    //2: it does work?
    Memory(void) : Buffer(NULL), size(0) {}

public:

    Memory(unsigned int size) : Buffer(NULL), size(size)
    {
        assert(size > 0);
        Buffer = new T[size];
    }

    Memory(const Memory &rhs) : Buffer(NULL), size(rhs.size)
    {
        Buffer = new T[size];
        assert(Buffer != NULL && rhs.Buffer != NULL);
        std::copy(rhs.Buffer, rhs.Buffer + size, Buffer);
    }

    bool operator==(const Memory &rhs) const
    {
        return std::equal(Buffer, Buffer + size, rhs.Buffer);
    }

    Memory &operator=(const Memory &rhs)
    {
        if (this != &rhs)
        {
            T *BufferTmp = new T[rhs.size]; //tmp var assures valid Buffer if new throws Exception
            delete[] Buffer;
            Buffer = BufferTmp;
            size = rhs.size;
            std::copy(rhs.Buffer, rhs.Buffer + size, Buffer);
            //faster then: memcpy((void *) Buffer, (const void *) rhs.Buffer, sizeof(T)*size);
        }
        return *this;
    }

    virtual ~Memory(void)
    {
        delete[] Buffer;
    }

    inline unsigned int Size(void)
    {
        return size;
    }

    inline const T *data(void) const
    {
        return Buffer;
    }

    inline T *data(void)
    {
        return Buffer;
    }

    //DEBUG
    inline bool IsDynamic(void) const
    {
        return true;
    }
};

#define MEMORY_WIDTH (TWidth == MemoryConst::Dynamic ? Width : TWidth)
#define MEMORY_HEIGHT (THeight == MemoryConst::Dynamic ? Height : THeight)
//merge with GenericBuffer from GameEngine
template <class T, int TWidth, int THeight> class GenericBuffer_dyn
{
protected:

    Memory<T, (TWidth == MemoryConst::Dynamic && THeight == MemoryConst::Dynamic) ? MemoryConst::Dynamic : TWidth * THeight> memory;
    COORDINATE Size;

public:

    GenericBuffer_dyn<T, TWidth, THeight>(void) : memory()
    {
        //if Buffer is dynamic, compiler error is thrown, which is correct
        //but no error is thrown?
        Size.X = TWidth;
        Size.Y = THeight;
    }

    GenericBuffer_dyn<T, TWidth, THeight>(unsigned int Width, unsigned int Height) : memory(Width*Height)
    {
        Size.X = MEMORY_WIDTH;
        Size.Y = MEMORY_HEIGHT;
        //assert(TWidth == Width && THeight == Height);
        /*std::tcout << "GenericBuffer_dyn::GenericBuffer_dyn(): Created <" << TWidth << "," << THeight << ">(" << Width << "," << Height << ")" << std::endl;
        std::tcout << "GenericBuffer_dyn::GenericBuffer_dyn(): Got " << GetWidth() << "x" << GetHeight() << " Matrix" << std::endl;*/
    }

    bool operator==(const GenericBuffer_dyn &rhs) const
    {
        return memory == rhs.memory;
    }

    unsigned int GetWidth(void) const
    {
        return Size.X;
    }

    unsigned int GetHeight(void) const
    {
        return Size.Y;
    }

    inline const T *Get(unsigned int x, unsigned int y) const
    {
        assert(x < Size.X && y < Size.Y);
        return memory.data() + y * Size.X + x;
    }

    inline T *Get(unsigned int x, unsigned int y)
    {
        assert(x < Size.X && y < Size.Y);
        return memory.data() + y * Size.X + x;
    }

    inline const T *PeriodicGet(int x, int y) const
    {
        //int % unsigned int == BAD!!!!!!!!!!!!!!!!!!!!!!
        return memory.data() + intmod(y, int(Size.Y)) * Size.X + intmod(x, int(Size.X));
    }

    inline T *PeriodicGet(int x, int y)
    {
        //int % unsigned int == BAD!!!!!!!!!!!!!!!!!!!!!!
        return memory.data() + intmod(y, int(Size.Y)) * Size.X + intmod(x, int(Size.X));
    }

    T operator[](unsigned int i) const
    {
        return memory.data()[i];
    }

    inline void Set(unsigned int x, unsigned int y, T Value)
    {
        assert(x < Size.X && y < Size.Y);
        memory.data()[y * Size.X + x] = Value;
    }

    inline bool IsDynamic(void)
    {
        return memory.IsDynamic();
    }

    void Zero(void)
    {
        std::memset(memory.data(), 0, Size.X * Size.Y * sizeof(T));
    }
};


#endif //LIBMEMORY_HEADER
