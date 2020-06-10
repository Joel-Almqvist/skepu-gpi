#pragma once
#ifndef SKEPU_STARPU_TENSOR3_PARTITION_HPP
#define SKEPU_STARPU_TENSOR3_PARTITION_HPP 1

#include <starpu_mpi.h>

#include <skepu3/cluster/common.hpp>
#include "../partition.hpp"

namespace skepu {
namespace util {

template<typename T>
class tensor3_partition : private partition_base<T>
{
	typedef partition_base<T> base;

	size_t m_size_i;
	size_t m_size_j;
	size_t m_size_k;
	size_t m_size_jk;
	size_t m_part_i;

public:
	typedef skepu::Index3D index_type;

	tensor3_partition() noexcept
	: base()
	{}

	tensor3_partition(size_t i, size_t j, size_t k) noexcept
	: base(), m_size_i(i), m_size_j(j), m_size_k(k), m_size_jk(j*k)
	{
		auto count = i * j * k;
		set_sizes();
		if(count)
			alloc_partitions();
	}

	tensor3_partition(tensor3_partition const & other) noexcept
	: tensor3_partition(other.m_size_i, other.m_size_j, other.m_size_k)
	{
		base::copy(other);
	}

	tensor3_partition(tensor3_partition && other) noexcept
	: base(std::move(other))
	{}

	~tensor3_partition() noexcept = default;

	auto
	operator=(tensor3_partition const & other) noexcept
	-> tensor3_partition &
	{
		this->~tensor3_partition();
		new(this) tensor3_partition(other);
		return *this;
	}

	auto
	operator=(tensor3_partition && other) noexcept
	-> tensor3_partition &
	{
		this->~tensor3_partition();
		new(this) tensor3_partition(std::move(other));
		return *this;
	}

	auto
	operator()(size_t i, size_t j, size_t k) noexcept
	-> T &
	{
		return base::operator()((i * m_size_jk) + (j * m_size_k) + k);
	}

	auto
	operator()(size_t i, size_t j, size_t k) const noexcept
	-> T const &
	{
		return base::operator()((i * m_size_jk) + (j * m_size_k) + k);
	}

	auto
	getParent() noexcept
	-> tensor3_partition &
	{
		return *this;
	}

	auto
	index(size_t pos) const noexcept
	-> index_type
	{
		index_type idx;
		idx.i = pos / m_size_jk;
		pos -= idx.i * m_size_jk;
		idx.j = pos / m_size_k;
		idx.k = pos - (idx.j * m_size_k);

		return idx;
	}

	using base::operator();
	using base::allgather;
	using base::block_count_from;
	using base::capacity;
	using base::data;
	using base::fill;
	using base::gather_to_root;
	using base::local_storage_handle;
	using base::handle_for;
	using base::partition;
	using base::scatter_from_root;
	using base::set;
	using base::size;

private:
	auto
	alloc_local_storage() noexcept
	-> void override
	{
		base::m_data = new T[base::m_size];
		starpu_block_data_register(
			&(base::m_data_handle),
			STARPU_MAIN_RAM,
			(uintptr_t)(base::m_data),
			m_size_k, m_size_k * m_size_j,
			m_size_k, m_size_j, m_size_i,
			sizeof(T));
		starpu_mpi_data_register(
			base::m_data_handle,
			cluster::mpi_tag(),
			STARPU_MAIN_RAM);
	}

	auto
	alloc_partitions() noexcept
	-> void override
	{
		if(base::m_size)
		{
			base::m_part_data = new T[base::m_part_size];
			for(size_t i(0); i < base::m_handles.size(); ++i)
			{
				auto ptr = (i == cluster::mpi_rank() ? base::m_part_data : 0);

				/* StarPU manages memory on nodes that does not own the data. */
				int home_node = -1;
				/* But it the rank owns the data, we manage the storage. */
				if(ptr)
					home_node = STARPU_MAIN_RAM;

				starpu_block_data_register(
					&base::m_handles[i],
					home_node,
					(uintptr_t)ptr,
					m_size_k, m_size_j * m_size_k,
					m_size_k, m_size_j, m_part_i,
					sizeof(T));
				starpu_mpi_data_register(
						base::m_handles[i],
						cluster::mpi_tag(),
						i);
			}
		}
	}

	auto
	get_ptr(starpu_data_handle_t & handle) noexcept
	-> T * override
	{
		return (T *)starpu_block_get_local_ptr(handle);
	}

	auto
	set_sizes() noexcept
	-> void
	{
		base::m_size = m_size_i * m_size_j * m_size_k;
		if(base::m_size)
		{
			auto ranks = skepu::cluster::mpi_size();
			m_part_i = m_size_i / ranks;
			if(m_size_i - (m_part_i * ranks))
				++m_part_i;
			base::m_part_size = m_part_i * m_size_j * m_size_k;
			base::m_capacity = base::m_part_size * ranks;
		}
	}
};

} // namespace util
} // namespace skepu

#endif // SKEPU_STARPU_TENSOR3_PARTITION_HPP
