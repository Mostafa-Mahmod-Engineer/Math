#pragma once
#include"Utilities.h"
#define fileName "Data Structures.h"


//class string______________________|
/*{
	int size = 0;
	char*buff = nullptr;
public:

};*/
//__________________________________|
template<class T, class H>
struct pair
{
	T first; H second;
	pair() {}
	pair(const T& f, const H& s) : first(f), second(s) {}
	pair(pair<T, H>&& other) { first = std::move(other.first); second = std::move(other.second); }
	pair<T, H>& operator=(pair<T, H>&& other)
	{
		pair<T, H> res;
		res.first = std::move(other.first);
		res.second = std::move(other.second);
		return res;
	}
	bool operator==(const pair<T, H>& other) { return this->first == other.first && this->second == other.second; }
	bool operator!=(const pair<T, H>& other) { return !(*this == other); }
	void print() { std::cout << " ( " << first << " , " << second << " ) "; }
	~pair() {}
};
template<class T, class H>
std::ostream& operator<<(std::ostream& C, const pair<T, H>& p)
{
	C << p.first << p.second;
	return C;
}
template<class T>
struct node
{
	T data;
	node<T> *next = nullptr, *prev = nullptr;
	int index = 0;
	node() {}
	node(T val) : data(val) {}
	node(T val, int i) : data(val), index(i) {}
	node(node<T>& other) { data = other.data; next = other.next; prev = other.prev; index = other.index; }
	node(node<T>* other) { data = other->data; next = other->next; prev = other->prev; index = other->index; }
	node<T>& operator=(node<T>& other)
	{
		this->data = other.data;
		this->next = other.next;
		this->prev = other.prev;
		this->index = other.index;
		return *this;
	}
	node<T>& operator=(node<T>* other)
	{
		this->data = other->data;
		this->next = other->next;
		this->prev = other->prev;
		this->index = other->index;
		return *this;
	}
	bool operator==(node<T>& other) { return this->data == other.data; }
	bool operator!=(node<T>& other) { return this->data != other.data; }
	node<T>& operator++() { if (this->next != nullptr)  *this = *(this->next); return *this; }
	node<T>& operator--() { if (this->prev != nullptr) *this = *(this->prev); return *this; }
	void print() { std::cout << this->data; }
	~node() {}
};

template<class T>
class list
{
	int size = 0;
	node<T> *head = nullptr, *tail = nullptr;
public:
	list() {}
	//______________________________________________________________________________
	template<class J>															//  |
	void assign(const J& j) { add(j); }											//  |
	template<class H, class ... J>												//  |
	void assign(const H& h, const J& ... j) { add(h); assign(j...); }			//  |
	template<class H, class ... J>												//  |
	list(const H& h, const J& ... j) { assign(h, j...); }						//  |
	//______________________________________________________________________________|
	list(const list<T>& other)
	{
		node<T>* temp = other.head;
		for (register int i = 0; i < other.size; i++)
		{
			this->add(temp->data);
			temp = temp->next;
		}
	}
	list(list<T>&& other)
	{
		this->size = other.size;
		this->head = other.head;
		this->tail = other.tail;
		other.size = 0;
		other.head = nullptr;
		other.tail = nullptr;
	}
	list(const list<T>* other)
	{
		node<T>* temp = other->head;
		for (register int i = 0; i < other->size; i++)
		{
			this->add(temp->data);
			temp = temp->next;
		}
	}
	list<T>& operator=(const list<T>& other)
	{
		this->~list();
		node<T>* temp = other.head;
		for (register int i = 0; i < other.size; i++)
		{
			this->add(temp->data);
			temp = temp->next;
		}
		return *this;
	}
	list<T>& operator=(list<T>&& other)
	{
		this->~list();
		this->size = other.size;
		this->head = other.head;
		this->tail = other.tail;
		other.size = 0;
		other.head = nullptr;
		other.tail = nullptr;
		return *this;//test the lvalue refrence copy
	}
	void insert_front(const T& val)
	{
		node<T>* newnode = new node<T>(val);
		if (size == 0)
		{
			head = tail = newnode;
			head->index = 0;
			tail->index = 0;
		}
		else
		{
			newnode->next = head;
			head->prev = newnode;
			head = newnode;
			node<T>* temp = head;
			for (register int i = 0; i < size; i++)
			{
				temp->index = i;
				temp = temp->next;
			}
		}
		size++;
	}
	void add(const T& val)
	{
		node<T>* newnode = new node<T>(val);
		if (size == 0)
		{
			head = tail = newnode;
			head->index = 0;
			tail->index = 0;
		}
		else if (size == 1)
		{
			tail = newnode;
			head->next = tail;
			tail->prev = head;
			tail->index = 1;
		}
		else
		{
			newnode->prev = tail;
			tail->next = newnode;
			tail = newnode;
			tail->index = size;
		}
		size++;
	}
	void insert(const T& val, int pos)
	{
		chk_rng(pos, size);
		if (pos == size)
			add(val);
		else
		{
			node<T>* temp = head;
			for (register int i = 0; i < pos; i++)
				temp = temp->next;
			node<T>* newnode = new node<T>(val);
			newnode->next = temp;
			newnode->prev = temp->prev;
			temp->prev = newnode;
			newnode->prev->next = newnode;
			size++;
			temp = head;
			for (register int i = 0; i < size; i++)
			{
				temp->index = i;
				temp = temp->next;
			}
		}
	}
	void remove_front()
	{
		head = head->next;
		delete head->prev;
		size--;
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
		{
			temp->index = i;
			temp = temp->next;
		}
	}
	void remove_last()
	{
		tail = tail->prev;
		delete tail->next;
		size--;
	}
	void remove(int pos)
	{
		chk_rng(pos, size);
		if (pos == size-1)
			remove_last();
		else if (pos == 0)
			remove_front();
		else
		{
			node<T>* temp;
			temp = head;
			for (register int i = 0; i < pos-1; i++)
				temp = temp->next;
			temp->next = temp->next->next;
			temp->next->prev = temp;
			size--;
			temp = head;
			for (register int i = 0; i < size; i++)
			{
				temp->index = i;
				temp = temp->next;
			}
		}
	}
	list<T> tie(const list<T>& other)
	{
		list<T> res1(this), res2(other);
		res1.tail->next = res2.head;
		res2.head->prev = res1.tail;
		res1.size += res2.size;
		node<T>* temp = res2.head, *Temp = res1.head;
		while (temp != NULL)
		{
			temp->index += Temp->index;
			temp = temp->next;
			Temp = Temp->next;
		}
		res1.tail = res2.tail;
		return res1;
	}
	T get(int pos)
	{
		if (pos > size)
			error("Position is Greater that the Size of the List");
		else if (pos < 0)
			error("Negative Index\nError Call from list::get()");
		else if (pos == 0)
			return head->data;
		else if (pos == size)
			return tail->data;
		if (pos < size / 2)
		{
			node<T>* temp = head;
			for (register int i = 0; i < pos; i++)
				temp = temp->next;
			return temp->next->data;
		}
		else
		{
			node<T>* temp = tail;
			for (register int i = 0; i < pos; i++)
				temp = temp->prev;
			return temp->prev->data;
		}
	}
	int find(const T& val)
	{
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
			if (temp->data == val)
				return i;
			else
				temp = temp->next;
		return -1;
	}
	void print()
	{
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
		{
			std::cout << std::endl << temp->index << ". ";
			std::cout << temp->data;
			temp = temp->next;
		}
		std::cout << std::endl;
	}
	void printback()
	{
		node<T>* temp = tail;
		for (register int i = size; i > 0; i--)
		{
			std::cout << std::endl << temp->index << ". ";
			std::cout << temp->data;
			temp = temp->prev;
		}
		std::cout << std::endl;
	}
	T Head() { return head->data; }
	T Tail() { return tail->data; }
	int Size() { return size; }
	node<T>* begin() { return head; }
	node<T>* end() { return tail; }

	list<T> tie(const list<T>& other) const
	{
		list<T> res1(this), res2(other);
		res1.tail->next = res2.head;
		res2.head->prev = res1.tail;
		res1.size += res2.size;
		node<T>* temp = res2.head, *Temp = res1.head;
		while (temp != NULL)
		{
			temp->index += Temp->index;
			temp = temp->next;
			Temp = Temp->next;
		}
		res1.tail = res2.tail;
		return res1;
	}
	T get(int pos) const
	{
		if (pos >= size)
		{
			printf("\nError\nPosition is Greater that the Size of the List");
			std::cin.ignore();
		}
		else if (pos == size - 1)
			return tail->data;
		if (pos < size / 2)
		{
			node<T>* temp = head;
			for (register int i = 0; i < pos; i++)
				temp = temp->next;
			return temp->data;
		}
		else
		{
			node<T>* temp = tail;
			for (register int i = 0; i < pos; i++)
				temp = temp->prev;
			return temp->data;
		}
	}
	int find(const T& val) const
	{
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
			if (temp->data == val)
				return i;
			else
				temp = temp->next;
		return -1;
	}
	void print() const
	{
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
		{
			std::cout << std::endl << temp->index << ". ";
			std::cout << temp->data;
			temp = temp->next;
		}
		std::cout << std::endl;
	}
	void printback() const
	{
		node<T>* temp = tail;
		for (register int i = size; i > 0; i--)
		{
			std::cout << std::endl << temp->index << ". ";
			std::cout << temp->data;
			temp = temp->prev;
		}
		std::cout << std::endl;
	}
	T Head() const { return head->data; }
	T Tail() const { return tail->data; }
	int Size() const { return size; }
	node<T>* begin() const { return head; }
	node<T>* end() const { return tail; }

	~list()
	{
		if (size)
		{
			if (size == 1)
				delete head;
			else if (size == 2)
			{
				delete head; delete tail;
			}
			else
			{
				node<T>* temp = head;
				do
				{
					temp = temp->next;
					delete temp->prev;
				} while (temp->next != nullptr);
				delete tail;
				temp = nullptr;
			}
		}
	}
};
template<class T>
class loop
{
	int size = 0;
	node<T> *head = nullptr, *tail = nullptr;
public:
	loop() {}
	loop(const loop<T>& other)
	{
		node<T>* temp = other.head;
		for (register int i = 0; i < other.size; i++)
		{
			this->add(temp->data);
			temp = temp->next;
		}
	}
	loop(loop<T>&& other)
	{
		this->size = other.size;
		this->head = other.head;
		this->tail = other.tail;
		other.size = 0;
		other.head = nullptr;
		other.tail = nullptr;
	}
	loop(const loop<T>* other)
	{
		node<T>* temp = other->head;
		for (register int i = 0; i < other->size; i++)
		{
			this->add(temp->data);
			temp = temp->next;
		}
	}
	loop<T>& operator=(const loop<T>& other)
	{
		this->~loop();
		node<T>* temp = other.head;
		for (register int i = 0; i < other.size; i++)
		{
			this->add(temp->data);
			temp = temp->next;
		}
		return *this;
	}
	loop<T>& operator=(loop<T>&& other)
	{
		this->~loop();
		this->size = other.size;
		this->head = other.head;
		this->tail = other.tail;
		other.size = 0;
		other.head = nullptr;
		other.tail = nullptr;
		return *this;
	}
	void insert_front(const T& val)
	{
		node<T>* newnode = new node<T>(val);
		if (size == 0)
		{
			head = tail = newnode;
			head->index = 0;
			tail->index = 0;
		}
		else if (size == 1)
		{
			newnode->next = head;
			head->prev = newnode;
			head = newnode;
			head->prev = tail;
			tail->next = head;
			tail->index = 1;
		}
		else
		{
			newnode->next = head;
			head->prev = newnode;
			head = newnode;
			node<T>* temp = head;
			for (register int i = 0; i < size; i++)
			{
				temp->index = i;
				temp = temp->next;
			}
		}
		size++;
	}
	void add(const T& val)
	{
		node<T>* newnode = new node<T>(val);
		if (size == 0)
		{
			head = tail = newnode;
			head->index = 0;
			tail->index = 0;
		}
		else if (size == 1)
		{
			tail = newnode;
			head->next = tail;
			tail->prev = head;
			head->prev = tail;
			tail->next = head;
			tail->index = 1;
		}
		else
		{
			newnode->prev = tail;
			tail->next = newnode;
			tail = newnode;
			tail->index = size;
		}
		size++;
	}
	void insert(const T& val, int pos)
	{
		if (pos > size)
		{
			printf("Error\nPosition Exceeds the Size of the loop");
			std::cin.ignore();
		}
		else if (pos == size)
			add(val);
		else if (pos == 0)
			insert_front(val);
		else
		{
			node<T>* temp = head;
			for (register int i = 0; i < pos; i++)
				temp = temp->next;
			node<T>* newnode = new node<T>(val);
			newnode->next = temp;
			newnode->prev = temp->prev;
			temp->prev = newnode;
			newnode->prev->next = newnode;
			size++;
			temp = temp->prev;
			for (register int i = pos; i < size; i++)
			{
				temp->index = i;
				temp = temp->next;
			}
		}
	}
	void remove_front()
	{
		head = head->next;
		delete head->prev;
		size--;
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
		{
			temp->index = i;
			temp = temp->next;
		}
	}
	void remove_last()
	{
		tail = tail->prev;
		delete tail->next;
		size--;
	}
	void remove(int pos)
	{
		if (pos > size)
		{
			printf("Error\nPosition Exceeds the Size of the loop");
			std::cin.ignore();
		}
		if (pos == size)
			remove_last();
		else
		{
			node<T>* temp;
			if (pos < size / 2)
			{
				temp = head;
				for (register int i = 0; i < pos; i++)
					temp = temp->next;
				temp->next = temp->next->next;
				temp->next->prev = temp;
			}
			else
			{
				temp = tail;
				for (register int i = size; i > pos; i--)
					temp = temp->prev;
				temp->prev = temp->prev->prev;
				temp->prev->next = temp;
			}
			size--;
			temp = temp->prev->prev;
			for (register int i = temp->index; i < size; i++)
			{
				temp->index = i;
				temp = temp->next;
			}
		}
	}
	loop<T> tie(const loop<T>& other)
	{
		loop<T> res1(this), res2(other);
		res1.tail->next = res2.head;
		res2.head->prev = res1.tail;
		res1.size += res2.size;
		node<T>* temp = res2.head, *Temp = res1.head;
		while (temp != NULL)
		{
			temp->index += Temp->index;
			temp = temp->next;
			Temp = Temp->next;
		}
		res1.tail = res2.tail;
		return res1;
	}
	T get(int pos)
	{
		if (pos > size)
		{
			printf("\nError\nPosition is Greater that the Size of the List");
			std::cin.ignore();
		}
		else if (pos == size)
			return tail->data;
		if (pos < size / 2)
		{
			node<T>* temp = head;
			for (register int i = 0; i < pos; i++)
				temp = temp->next;
			return temp->next->data;
		}
		else
		{
			node<T>* temp = tail;
			for (register int i = 0; i < pos; i++)
				temp = temp->prev;
			return temp->prev->data;
		}
	}
	int find(const T& val)
	{
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
			if (temp->data == val)
				return i;
			else
				temp = temp->next;
		return -1;
	}
	void print()
	{
		node<T>* temp = head;
		for (register int i = 0; i < size; i++)
		{
			std::cout << std::endl << temp->index << ". ";
			std::cout << temp->data;
			temp = temp->next;
		}
		std::cout << std::endl;
	}
	void printback()
	{
		node<T>* temp = tail;
		for (register int i = size; i > 0; i--)
		{
			std::cout << std::endl << temp->index << ". ";
			std::cout << temp->data;
			temp = temp->prev;
		}
		std::cout << std::endl;
	}
	T Head() { return head->data; }
	T Tail() { return tail->data; }
	node<T>* H_iter() { return head; }
	node<T>* T_iter() { return tail; }
	~loop()
	{
		if (size)
		{
			node<T>* temp = head;
			do
			{
				temp = temp->next;
				delete temp->prev;
			} while (temp->next != nullptr);
			delete tail;
			temp = nullptr;
		}
	}
};
template<class T>
class array
{
	T* arr = nullptr; int size, last; bool sorted = false;
public:
	array(int s) : size(s) { arr = new T[size]; last = -1; }
	template<class G>
	void assign(const G& g) { arr[size - 1] = g; }
	template<class H, class ... G>
	void assign(const H& h, const G& ... g)
	{
		arr[size - sizeof...(g) - 1] = h;
		assign(g...);
	}
	template<class H, class ... G>
	array(const H& h, const G& ... g)
	{
		size = sizeof...(G) + 1;
		arr = new T[size];
		assign(h, g...);
		last = size - 1;
	}
	array(const array<T>& other)
	{
		size = other.size;
		arr = new T[size];
		for (register int i = 0; i < size; i++)
			arr[i] = other.arr[i];
		last = other.last;
		sorted = other.sorted;
	}
	array(array<T>&& other)
	{
		size = other.size;
		arr = other.arr;
		last = other.last;
		sorted = other.sorted;
		other.arr = nullptr;
	}
	array(const array<T>* other)
	{
		size = other->size;
		arr = new T[size];
		for (register int i = 0; i < size; i++)
			arr[i] = other->arr[i];
		last = other->last;
		sorted = other->sorted;
	}
	array<T>& operator=(const array<T>& other)
	{
		if (arr == other.arr)
			return *this;
		this->~array();
		size = other.size;
		arr = new T[size];
		for (register int i = 0; i < size; i++)
			arr[i] = other.arr[i];
		last = other.last;
		sorted = other.sorted;
		return *this;
	}
	array<T>& operator=(array<T>&& other)
	{
		size = other.size;
		arr = other.arr;
		last = other.last;
		sorted = other.sorted;
		other.arr = nullptr;
		return *this;

	}
	void resize(int s)
	{
		size = s;
		T* newarr = new T[size];
		for (register int i = 0; i < last; i++)
			newarr[i] = arr[i];
		delete[] arr;
		arr = newarr;
		newarr = nullptr;
	}
	void resize()
	{
		size++;
		T* newarr = new T[size];
		for (register int i = 0; i < last; i++)
			newarr[i] = arr[i];
		delete[] arr;
		arr = newarr;
		newarr = nullptr;
	}
	void reverse()
	{
		T temp;
		for (register int i = 0; i < last / 2; i++)
		{
			temp = arr[i];
			arr[i] = arr[last - i];
			arr[last - i] = temp;
		}
	}
	/* impliment temp to avoid copying ====>>*/void _merge(int start, int mid, int end)
	{
		list<T> temp;
		int i = start, j = mid + 1;
		while (i <= mid && j <= end)
			if (arr[i] > arr[j])
			{
				temp.add(arr[j]);
				j++;
			}
			else
			{
				temp.add(arr[i]);
				i++;
			}
		while (i <= mid)
		{
			temp.add(arr[i]);
			i++;
		}
		while (j <= end)
		{
			temp.add(arr[j]);
			j++;
		}
		node<T>* n = temp.begin();
		for (i = start; i < temp.Size() + start; i++, n = n->next)
			this->arr[i] = n->data;
	}
	void sort(int start = 0, int end = -1)
	{
		sorted = true;
		if (end == -1)
			end = last;
		if (start < end)
		{
			int mid = start / 2.0 + end / 2.0;
			sort(start, mid);
			sort(mid + 1, end);
			_merge(start, mid, end);
		}
		return;
	}
	void push(const T& val)
	{
		if (last == size - 1)
		{
			printf("\nError\nArray is Full");
			std::cin.ignore();
		}
		last++;
		arr[last] = val;
	}
	void push_back(const T& val)
	{
		if (last == size - 1)
			resize();
		last++;
		arr[last] = val;
	}
	T pop_back()
	{
		if (last == 0)
		{
			printf("\nError\nArray is Emtpy");
			std::cin.ignore();
		}
		T res = arr[last];
		delete arr[last];
		last--;
		return res;

	}
	T peek() { return arr[last]; }
	T& operator[](int i) { return arr[i]; }
	T* begin() { return arr; }
	T* end() { return arr + sizeof(T)*last; }
	int Size() { return size; }
	int length() { return last; }
	int find(const T& val)
	{
		for (register int i = 0; i < last; i++)
			if (arr[i] == val)
				return i;
		return -1;
	}
	T accumulate()
	{
		T res = arr[0];
		for (register int i = 1; i < last; i++)
			res += arr[i];
		return res;
	}
	bool is_sorted()
	{
		if (sorted)
			return true;
		int count = 0;
		while (arr[count] == arr[count + 1])
			count++;
		if (arr[count] < arr[count + 1])
			for (register int i = count + 2; i < last; i++)
				if (arr[i - 1] > arr[i])
					return false;
				else
					continue;
		else
			for (register int i = count + 2; i < last; i++)
				if (arr[i - 1] < arr[i])
					return false;
		sorted = true;
		return true;
	}
	int binary_search(const T& val)
	{
		if (is_sorted())
		{
			int start = 0, end = last, mid;
			while (start != end)
			{
				mid = start / 2.0 + end / 2.0;
				if (val == arr[mid])
					return mid;
				else if (val > arr[mid])
					start = mid;
				else
					end = mid;
			}
			return -1;
		}
		printf("\nError\nArray not Sorted"); while (1);
	}
	void print()
	{
		printf("[");
		for (register int i = 0; i < last; i++)
			std::cout << arr[i] << " , ";
		std::cout << arr[last] << "]\n";
	}

	void resize(int s) const
	{
		size = s;
		T* newarr = new T[size];
		for (register int i = 0; i < last; i++)
			newarr[i] = arr[i];
		delete[] arr;
		arr = newarr;
		newarr = nullptr;
	}
	void resize() const
	{
		size++;
		T* newarr = new T[size];
		for (register int i = 0; i < last; i++)
			newarr[i] = arr[i];
		delete[] arr;
		arr = newarr;
		newarr = nullptr;
	}
	void reverse() const
	{
		T temp;
		for (register int i = 0; i < last / 2; i++)
		{
			temp = arr[i];
			arr[i] = arr[last - i];
			arr[last - i] = temp;
		}
	}
	T peek() const { return arr[last]; }
	T& operator[](int i) const { return arr[i]; }
	T* begin() const { return arr; }
	T* end() const { return arr + sizeof(T)*last; }
	int Size() const { return size; }
	int length() const { return last; }
	int find(const T& val) const
	{
		for (register int i = 0; i < last; i++)
			if (arr[i] == val)
				return i;
		return -1;
	}
	T accumulate() const
	{
		T res = arr[0];
		for (register int i = 1; i < last; i++)
			res += arr[i];
		return res;
	}
	bool is_sorted() const
	{
		if (sorted)
			return true;
		int count = 0;
		while (arr[count] == arr[count + 1])
			count++;
		if (arr[count] < arr[count + 1])
			for (register int i = count + 1; i < last; i++)
				if (arr[i - 1] > arr[i])
					return false;
				else
					continue;
		else
			for (register int i = count + 1; i < last; i++)
				if (arr[i - 1] < arr[i])
					return false;
		sorted = true;
		return true;
	}
	int binary_search(const T& val) const
	{
		if (is_sorted())
		{
			int start = 0, end = last, mid;
			while (start != end)
			{
				mid = start / 2.0 + end / 2.0;
				if (val == arr[mid])
					return mid;
				else if (val > arr[mid])
					start = mid;
				else
					end = mid;
			}
			return -1;
		}
		printf("\nError\nArray not Sorted"); while (1);
	}
	void print() const
	{
		printf("[");
		for (register int i = 0; i < last; i++)
			std::cout << arr[i] << " , ";
		std::cout << arr[last] << "]\n";
	}

	~array() { if (arr != nullptr) delete[] arr; }
};
template<class T>
class stack
{
	int size = 0, top;
	T* stc;
public:
	stack(int s) : size(s) { stc = new T[s]; }
	stack(const stack<T>& other)
	{
		size = other.size;
		top = other.top;
		stc = new T[size];
		for (register int i = 0; i <= top; i++)
			stc[i] = other.stc[i];
	}
	stack(stack<T>&& other)
	{
		size = other.size;
		top = other.top;
		stc = other.stc;
		other.stc = nullptr;
	}
	stack(const stack<T>* other)
	{
		size = other->size;
		top = other->top;
		stc = new T[size];
		for (register int i = 0; i <= top; i++)
			stc[i] = other->stc[i];
	}
	stack<T>& operator=(const stack<T>& other)
	{
		size = other.size;
		top = other.top;
		stc = new T[size];
		for (register int i = 0; i <= top; i++)
			stc[i] = other.stc[i];
		return *this;
	}
	stack<T>& operator=(stack<T>&& other)
	{
		size = other.size;
		top = other.top;
		stc = other.stc;
		other.stc = nullptr;
		return *this;
	}
	void push(const T& val)
	{
		if (full())
		{
			printf("\nError\nStack is Full");
			std::cin.ignore();
		}
		top++;
		stc[top] = val;
	}
	T pop()
	{
		if (empty())
		{
			printf("\nError\nStack is Empty");
			std::cin.ignore();
		}
		T res = stc[top];
		top--;
		return res;
	}
	T peek()
	{
		if (empty())
		{
			printf("\nError\nStack is Empty");
			std::cin.ignore();
		}
		return stc[top];
	}
	bool full() { return top == size; }
	bool empty() { return top == 0; }
	void print()
	{
		for (register int i = size; i > 0; i--)
			std::cout << std::endl << i << ". " << stc[i];
	}
	~stack() { delete[] stc; }
};
template<class T>
class queue
{
	int size, back = -1;
	T* q;
public:
	queue(int s) : size(s) { q = new T[size]; }
	queue(const queue<T>& other)
	{
		size = other.size;
		back = other.back;
		q = new T[size];
		for (register int i = 0; i < back; i++)
			q[i] = other.q[i];
	}
	queue(queue<T>&& other)
	{
		size = other.size;
		back = other.back;
		q = other.q;
		other.q = nullptr;
	}
	queue(const queue<T>* other)
	{
		size = other->size;
		back = other->back;
		q = new T[size];
		for (register int i = 0; i < back; i++)
			q[i] = other->q[i];
	}
	queue<T>& operator=(const queue<T>& other)
	{
		if (this->q == other.q)
			return *this;
		this->~queue();
		this->size = other.size;
		this->back = other.back;
		this->q = new T[size];
		for (register int i = 0; i < back; i++)
			this->q[i] = other.q[i];
		return *this;
	}
	queue<T>& operator=(queue<T>&& other)
	{
		this->~queue();
		this->size = other.size;
		this->back = other.back;
		this->q = other.q;
		other.q = nullptr;
		return *this;
	}
	void enqueue(const T& val)
	{
		if (full())
		{
			printf("\nError\nQueue is Full");
			std::cin.ignore();
		}
		back++;
		q[back] = val;
	}
	T dequeue()
	{
		if (empty())
		{
			printf("\nError\nQueue is Empty");
			std::cin.ignore();
		}
		T res = q[0];
		for (register int i = 0; i < back; i++)
			q[i] = q[i + 1];
		back--;
		return res;
	}
	T peek()
	{
		if (empty())
		{
			printf("\nError\nQueue is Empty");
			std::cin.ignore();
		}
		return q[0];
	}
	bool full() { return back == size; }
	bool empty() { return back == 0; }
	void print()
	{
		for (int i = 0; i <= back; i++)
			std::cout << std::endl << i << ". " << q[i];
	}
	~queue() { delete[] q; }
};
template<class T>
class ostack
{
	list<T> stc;
public:
	ostack() {}
	void push(const T& val) { stc.add(val); }
	T pop()
	{
		if (stc.Size() == 0)
		{
			printf("\nError\nStack is Empty");
			std::cin.ignore();
		}
		T res = stc.Tail();
		stc.remove_last();
		return res;
	}
	T peek()
	{
		if (stc.Size() == 0)
		{
			printf("\nError\nStack is Empty");
			std::cin.ignore();
		}
		return stc.Tail();
	}
	void print() { stc.printback(); }
	~ostack() {}
};
template<class T>
class oqueue
{
	list<T> q;
public:
	oqueue() {}
	void enqueue(const T& val) { q.add(val); }
	T dequeue()
	{
		if (q.Size() == 0)
		{
			printf("\nError\nQueue is Empty");
			std::cin.ignore();
		}
		T res = q.Head();
		q.remove_front();
		return res;
	}
	T peek()
	{
		if (q.Size() == 0)
		{
			printf("\nError\nStack is Empty");
			std::cin.ignore();
		}
		return q.Head();
	}
	void print() { q.print(); }
	~oqueue() {}
};



#undef fileName