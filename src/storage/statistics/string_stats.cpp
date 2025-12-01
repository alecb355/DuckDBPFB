#include "duckdb/storage/statistics/string_stats.hpp"

#include "duckdb/common/serializer/deserializer.hpp"
#include "duckdb/common/serializer/serializer.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/main/error_manager.hpp"
#include "duckdb/storage/statistics/base_statistics.hpp"
#include "utf8proc_wrapper.hpp"
#include "duckdb/common/types/blob.hpp"

namespace duckdb {

BaseStatistics StringStats::CreateUnknown(LogicalType type) {
	BaseStatistics result(std::move(type));
	result.InitializeUnknown();
	auto &string_data = StringStats::GetDataUnsafe(result);
	StringStats::Init_PBF(string_data);
	for (idx_t i = 0; i < StringStatsData::MAX_STRING_MINMAX_SIZE; i++) {
		string_data.min[i] = 0;
		string_data.max[i] = 0xFF;
	}
	string_data.max_string_length = 0;
	string_data.has_max_string_length = false;
	string_data.has_unicode = true;
	return result;
}

BaseStatistics StringStats::CreateEmpty(LogicalType type) {
	BaseStatistics result(std::move(type));
	result.InitializeEmpty();
	auto &string_data = StringStats::GetDataUnsafe(result);
	StringStats::Init_PBF(string_data);
	for (idx_t i = 0; i < StringStatsData::MAX_STRING_MINMAX_SIZE; i++) {
		string_data.min[i] = 0xFF;
		string_data.max[i] = 0;
	}
	string_data.max_string_length = 0;
	string_data.has_max_string_length = true;
	string_data.has_unicode = false;
	return result;
}

StringStatsData &StringStats::GetDataUnsafe(BaseStatistics &stats) {
	D_ASSERT(stats.GetStatsType() == StatisticsType::STRING_STATS);
	return stats.stats_union.string_data;
}

const StringStatsData &StringStats::GetDataUnsafe(const BaseStatistics &stats) {
	D_ASSERT(stats.GetStatsType() == StatisticsType::STRING_STATS);
	return stats.stats_union.string_data;
}

bool StringStats::HasMaxStringLength(const BaseStatistics &stats) {
	if (stats.GetType().id() == LogicalTypeId::SQLNULL) {
		return false;
	}
	return StringStats::GetDataUnsafe(stats).has_max_string_length;
}

uint32_t StringStats::MaxStringLength(const BaseStatistics &stats) {
	if (!HasMaxStringLength(stats)) {
		throw InternalException("MaxStringLength called on statistics that does not have a max string length");
	}
	return StringStats::GetDataUnsafe(stats).max_string_length;
}

bool StringStats::CanContainUnicode(const BaseStatistics &stats) {
	if (stats.GetType().id() == LogicalTypeId::SQLNULL) {
		return true;
	}
	return StringStats::GetDataUnsafe(stats).has_unicode;
}

string GetStringMinMaxValue(const data_t data[]) {
	idx_t len;
	for (len = 0; len < StringStatsData::MAX_STRING_MINMAX_SIZE; len++) {
		if (!data[len]) {
			break;
		}
	}
	return string(const_char_ptr_cast(data), len);
}

string StringStats::Min(const BaseStatistics &stats) {
	return GetStringMinMaxValue(StringStats::GetDataUnsafe(stats).min);
}

string StringStats::Max(const BaseStatistics &stats) {
	return GetStringMinMaxValue(StringStats::GetDataUnsafe(stats).max);
}

void StringStats::ResetMaxStringLength(BaseStatistics &stats) {
	StringStats::GetDataUnsafe(stats).has_max_string_length = false;
}

void StringStats::SetMaxStringLength(BaseStatistics &stats, uint32_t length) {
	auto &data = StringStats::GetDataUnsafe(stats);
	data.has_max_string_length = true;
	data.max_string_length = length;
}

void StringStats::SetContainsUnicode(BaseStatistics &stats) {
	StringStats::GetDataUnsafe(stats).has_unicode = true;
}

void StringStats::Serialize(const BaseStatistics &stats, Serializer &serializer) {
	auto &string_data = StringStats::GetDataUnsafe(stats);
	serializer.WriteProperty(200, "min", string_data.min, StringStatsData::MAX_STRING_MINMAX_SIZE);
	serializer.WriteProperty(201, "max", string_data.max, StringStatsData::MAX_STRING_MINMAX_SIZE);
	serializer.WriteProperty(202, "has_unicode", string_data.has_unicode);
	serializer.WriteProperty(203, "has_max_string_length", string_data.has_max_string_length);
	serializer.WriteProperty(204, "max_string_length", string_data.max_string_length);
}

void StringStats::Deserialize(Deserializer &deserializer, BaseStatistics &base) {
	auto &string_data = StringStats::GetDataUnsafe(base);
	deserializer.ReadProperty(200, "min", string_data.min, StringStatsData::MAX_STRING_MINMAX_SIZE);
	deserializer.ReadProperty(201, "max", string_data.max, StringStatsData::MAX_STRING_MINMAX_SIZE);
	deserializer.ReadProperty(202, "has_unicode", string_data.has_unicode);
	deserializer.ReadProperty(203, "has_max_string_length", string_data.has_max_string_length);
	deserializer.ReadProperty(204, "max_string_length", string_data.max_string_length);
}

static int StringValueComparison(const_data_ptr_t data, idx_t len, const_data_ptr_t comparison) {
	for (idx_t i = 0; i < len; i++) {
		if (data[i] < comparison[i]) {
			return -1;
		} else if (data[i] > comparison[i]) {
			return 1;
		}
	}
	return 0;
}

static void ConstructValue(const_data_ptr_t data, idx_t size, data_t target[]) {
	idx_t value_size = size > StringStatsData::MAX_STRING_MINMAX_SIZE ? StringStatsData::MAX_STRING_MINMAX_SIZE : size;
	memcpy(target, data, value_size);
	for (idx_t i = value_size; i < StringStatsData::MAX_STRING_MINMAX_SIZE; i++) {
		target[i] = '\0';
	}
}

void StringStats::Update(BaseStatistics &stats, const string_t &value) {
	auto data = const_data_ptr_cast(value.GetData());
	auto size = value.GetSize();

	//! we can only fit 8 bytes, so we might need to trim our string
	// construct the value
	data_t target[StringStatsData::MAX_STRING_MINMAX_SIZE];
	ConstructValue(data, size, target);


	auto &string_data = StringStats::GetDataUnsafe(stats);
	// Prefix Bloom Filter Update Logic    
    // 1. Ensure PBF is initialized 
    if (!string_data.has_pbf) {
        Init_PBF(string_data);
    }
	// 2. Compute and set bits for prefixes at 1, 2, 4, 8 bytes
    const uint32_t prefix_lengths[] = {1, 2, 4, 8};

    for (int i = 0; i < StringStatsData::NUM_PREFIXES; i++) {
        uint32_t len = prefix_lengths[i];

        // We can only process if the string is long enough 
        if (size >= len) {
            // Hash the prefix using DuckDB's standard Hash
            hash_t hash_val = Hash(reinterpret_cast<const char*>(data), len); // cast unsigned char* and char*

            // Map hash to bit index
            uint32_t bit_index = hash_val % StringStatsData::NUM_BITS;
            uint32_t byte_idx = bit_index / 8;
            uint32_t bit_offset = bit_index % 8;

            // Set the bit 
            string_data.prefixes[i].bits[byte_idx] |= (1 << bit_offset);
        }
    }
	// update the min and max
	if (StringValueComparison(target, StringStatsData::MAX_STRING_MINMAX_SIZE, string_data.min) < 0) {
		memcpy(string_data.min, target, StringStatsData::MAX_STRING_MINMAX_SIZE);
	}
	if (StringValueComparison(target, StringStatsData::MAX_STRING_MINMAX_SIZE, string_data.max) > 0) {
		memcpy(string_data.max, target, StringStatsData::MAX_STRING_MINMAX_SIZE);
	}
	if (size > string_data.max_string_length) {
		string_data.max_string_length = UnsafeNumericCast<uint32_t>(size);
	}
	if (stats.GetType().id() == LogicalTypeId::VARCHAR && !string_data.has_unicode) {
		auto unicode = Utf8Proc::Analyze(const_char_ptr_cast(data), size);
		if (unicode == UnicodeType::UTF8) {
			string_data.has_unicode = true;
		} else if (unicode == UnicodeType::INVALID) {
			throw ErrorManager::InvalidUnicodeError(string(const_char_ptr_cast(data), size),
			                                        "segment statistics update");
		}
	}
}

void StringStats::SetMin(BaseStatistics &stats, const string_t &value) {
	ConstructValue(const_data_ptr_cast(value.GetData()), value.GetSize(), GetDataUnsafe(stats).min);
}

void StringStats::SetMax(BaseStatistics &stats, const string_t &value) {
	ConstructValue(const_data_ptr_cast(value.GetData()), value.GetSize(), GetDataUnsafe(stats).max);
}

void StringStats::Merge(BaseStatistics &stats, const BaseStatistics &other) {
	if (other.GetType().id() == LogicalTypeId::VALIDITY) {
		return;
	}
	if (other.GetType().id() == LogicalTypeId::SQLNULL) {
		return;
	}
	auto &string_data = StringStats::GetDataUnsafe(stats);
	auto &other_data = StringStats::GetDataUnsafe(other);

	// NEW: Prefix Bloom Filter Merge Logic
    // If the other segment has PBF data, we must merge it 
    if (other_data.has_pbf) {
        // If the destination doesn't have PBF yet, initialize it now
        if (!string_data.has_pbf) {
            Init_PBF(string_data);
        }

        // Merge all four prefix levels (1, 2, 4, 8 bytes)
        for (int i = 0; i < StringStatsData::NUM_PREFIXES; i++) {
            // We iterate over the bytes of the bitset and OR them together 
            // NUM_BITS is in bits, so we divide by 8 for bytes
            for (int j = 0; j < StringStatsData::NUM_BITS / 8; j++) {
                string_data.prefixes[i].bits[j] |= other_data.prefixes[i].bits[j];
            }
        }
    }
	if (StringValueComparison(other_data.min, StringStatsData::MAX_STRING_MINMAX_SIZE, string_data.min) < 0) {
		memcpy(string_data.min, other_data.min, StringStatsData::MAX_STRING_MINMAX_SIZE);
	}
	if (StringValueComparison(other_data.max, StringStatsData::MAX_STRING_MINMAX_SIZE, string_data.max) > 0) {
		memcpy(string_data.max, other_data.max, StringStatsData::MAX_STRING_MINMAX_SIZE);
	}
	string_data.has_unicode = string_data.has_unicode || other_data.has_unicode;
	string_data.has_max_string_length = string_data.has_max_string_length && other_data.has_max_string_length;
	string_data.max_string_length = MaxValue<uint32_t>(string_data.max_string_length, other_data.max_string_length);
}

FilterPropagateResult StringStats::CheckZonemap(const BaseStatistics &stats,
                                                ExpressionType comparison_type,
                                                array_ptr<const Value> constants) {
	auto &string_data = StringStats::GetDataUnsafe(stats);

	std::cout << "[PBF] Enter CheckZonemap\n";

	// 1) First run PBF for all constants
	bool any_pbf_used = false;        // at least one constant actually went through PBF
	bool any_pbf_no_pruning = false;  // at least one constant returned NO_PRUNING_POSSIBLE from PBF

	// These constants are "not rejected by PBF" and will be passed to the min/max zonemap
	std::vector<const Value *> survivors;

	for (auto &constant_value : constants) {
		D_ASSERT(constant_value.type() == stats.GetType());
		D_ASSERT(!constant_value.IsNull());
		auto &constant = StringValue::Get(constant_value);

		std::cout << "[PBF]   constant = \"" << constant << "\"\n";

		auto prefix_query = StringStats::GetPrefixCandidates(comparison_type, constant);
		std::cout << "[PBF]   GetPrefixCandidates returned "
		          << prefix_query.prefixes.size() << " prefix(es)\n";

		// No usable prefix: PBF is effectively a no-op, so we go straight to min/max
		if (prefix_query.prefixes.empty()) {
			any_pbf_no_pruning = true;
			survivors.push_back(&constant_value);
			std::cout << "[PBF]   no usable prefixes -> fall back to min/max zonemap\n";
			continue;
		}

		any_pbf_used = true;
		auto prefix_res = StringStats::CheckPBF(stats, prefix_query);

		std::cout << "[PBF]   CheckPBF result = "
		          << (prefix_res == FilterPropagateResult::FILTER_ALWAYS_FALSE
		              ? "FILTER_ALWAYS_FALSE"
		              : "NO_PRUNING_POSSIBLE")
		          << "\n";

		if (prefix_res == FilterPropagateResult::FILTER_ALWAYS_FALSE) {
			// This constant is ruled out by PBF: this segment cannot match this constant
			std::cout << "[PBF]   this constant ruled out by PBF -> skip min/max for it\n";
			// Do not add to survivors
			continue;
		} else {
			// PBF cannot prune this constant -> forward it to min/max
			any_pbf_no_pruning = true;
			std::cout << "[PBF]   PBF says NO_PRUNING_POSSIBLE -> keep for min/max\n";
			survivors.push_back(&constant_value);
		}
	}

	// If at least one constant used PBF and all those constants were rejected by PBF,
	// then this segment can never satisfy the predicate and can be fully pruned.
	if (any_pbf_used && !any_pbf_no_pruning) {
		std::cout << "[PBF]   all prefix-eligible constants ruled out by PBF -> "
		          << "FILTER_ALWAYS_FALSE for this segment\n";
		return FilterPropagateResult::FILTER_ALWAYS_FALSE;
	}

	// 2) Run the original min/max zonemap logic only on constants that PBF did not reject
	for (auto *val_ptr : survivors) {
		const auto &constant_value = *val_ptr;
		auto &constant = StringValue::Get(constant_value);

		std::cout << "[ZONEMAP]   running min/max for constant = \"" << constant << "\"\n";

		auto prune_result = CheckZonemap(string_data.min, StringStatsData::MAX_STRING_MINMAX_SIZE,
		                                 string_data.max, StringStatsData::MAX_STRING_MINMAX_SIZE,
		                                 comparison_type, constant);

		if (prune_result == FilterPropagateResult::NO_PRUNING_POSSIBLE) {
			std::cout << "[ZONEMAP]   min/max: NO_PRUNING_POSSIBLE\n";
			return FilterPropagateResult::NO_PRUNING_POSSIBLE;
		} else if (prune_result == FilterPropagateResult::FILTER_ALWAYS_TRUE) {
			std::cout << "[ZONEMAP]   min/max: FILTER_ALWAYS_TRUE\n";
			return FilterPropagateResult::FILTER_ALWAYS_TRUE;
		} else {
			std::cout << "[ZONEMAP]   min/max: FILTER_ALWAYS_FALSE for this constant\n";
		}
	}

	// If we reach this point:
	//   - either all constants were rejected by PBF, or
	//   - PBF was effectively a no-op and all surviving constants were rejected by min/max
	std::cout << "[ZONEMAP]   final: FILTER_ALWAYS_FALSE\n";
	return FilterPropagateResult::FILTER_ALWAYS_FALSE;
}


FilterPropagateResult StringStats::CheckZonemap(const_data_ptr_t min_data, idx_t min_len, const_data_ptr_t max_data,
                                                idx_t max_len, ExpressionType comparison_type, const string &constant) {
	auto data = const_data_ptr_cast(constant.c_str());
	idx_t size = constant.size();

	int min_comp = StringValueComparison(data, MinValue(min_len, size), min_data);
	int max_comp = StringValueComparison(data, MinValue(max_len, size), max_data);
	switch (comparison_type) {
	case ExpressionType::COMPARE_EQUAL:
	case ExpressionType::COMPARE_NOT_DISTINCT_FROM:
		if (min_comp >= 0 && max_comp <= 0) {
			return FilterPropagateResult::NO_PRUNING_POSSIBLE;
		} else {
			return FilterPropagateResult::FILTER_ALWAYS_FALSE;
		}
	case ExpressionType::COMPARE_NOTEQUAL:
	case ExpressionType::COMPARE_DISTINCT_FROM:
		if (min_comp < 0 || max_comp > 0) {
			return FilterPropagateResult::FILTER_ALWAYS_TRUE;
		}
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	case ExpressionType::COMPARE_GREATERTHANOREQUALTO:
	case ExpressionType::COMPARE_GREATERTHAN:
		if (max_comp <= 0) {
			return FilterPropagateResult::NO_PRUNING_POSSIBLE;
		} else {
			return FilterPropagateResult::FILTER_ALWAYS_FALSE;
		}
	case ExpressionType::COMPARE_LESSTHAN:
	case ExpressionType::COMPARE_LESSTHANOREQUALTO:
		if (min_comp >= 0) {
			return FilterPropagateResult::NO_PRUNING_POSSIBLE;
		} else {
			return FilterPropagateResult::FILTER_ALWAYS_FALSE;
		}
	default:
		throw InternalException("Expression type not implemented for string statistics zone map");
	}
}

FilterPropagateResult StringStats::CheckPBF(const BaseStatistics &stats, const PrefixQuery &query) {
	auto &string_data = StringStats::GetDataUnsafe(stats);
	std::cout << "[PBF] CheckPBF called: has_pbf=" << (string_data.has_pbf ? 1 : 0)
	          << ", num_prefixes=" << query.prefixes.size() << std::endl;

	if (!string_data.has_pbf || query.prefixes.empty()) {
		std::cout << "[PBF]   no PBF or empty prefix list -> NO_PRUNING_POSSIBLE\n";
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}

	std::cout << "[PBF]   prefixes: ";
	for (auto &c : query.prefixes) {
		std::cout << "\"" << c << "\" ";
	}
	std::cout << std::endl;

	// Pick the "strongest" candidate: the prefix with the longest length
	const std::string *best_candidate = &query.prefixes[0];
	for (auto &c : query.prefixes) {
		if (c.size() > best_candidate->size()) {
			best_candidate = &c;
		}
	}

	std::cout << "[PBF]   using strongest prefix: \"" << *best_candidate << "\"\n";

	auto data = const_data_ptr_cast(best_candidate->c_str());
	auto size = best_candidate->size();
	if (size == 0) {
		// Empty prefix cannot be used for pruning
		std::cout << "[PBF]   strongest prefix is empty -> NO_PRUNING_POSSIBLE\n";
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}

	bool any_level_checked = false;

	for (idx_t i = 0; i < StringStatsData::NUM_PREFIXES; i++) {
		auto prefix_length = static_cast<uint32_t>(string_data.prefixes[i].level);

		// Only check the level whose prefix length matches the candidate length
		if (prefix_length == 0 || prefix_length != size) {
			continue;
		}
		any_level_checked = true;

		auto hash_val = Hash(reinterpret_cast<const char *>(data), prefix_length);
		auto bit_index = hash_val % StringStatsData::NUM_BITS;
		auto byte_idx = bit_index / 8;
		auto bit_mask = static_cast<uint8_t>(1 << (bit_index % 8));

		std::cout << "[PBF]   level=" << prefix_length
		          << " bit_index=" << bit_index
		          << " byte_idx=" << byte_idx
		          << " bit_mask=" << (int)bit_mask << std::endl;

		if (string_data.prefixes[i].bits[byte_idx] & bit_mask) {
			// This segment may contain strings with this prefix -> cannot prune
			std::cout << "[PBF]   bit is set -> NO_PRUNING_POSSIBLE\n";
			return FilterPropagateResult::NO_PRUNING_POSSIBLE;
		} else {
			std::cout << "[PBF]   bit is NOT set at this level\n";
		}
	}

	if (!any_level_checked) {
		std::cout << "[PBF]   no matching prefix length level -> NO_PRUNING_POSSIBLE\n";
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}

	// None of the levels for this prefix length have the bit set -> safe to prune this segment
	std::cout << "[PBF]   no bits set for this prefix -> FILTER_ALWAYS_FALSE\n";
	return FilterPropagateResult::FILTER_ALWAYS_FALSE;
}

static uint32_t GetValidMinMaxSubstring(const_data_ptr_t data) {
	for (uint32_t i = 0; i < StringStatsData::MAX_STRING_MINMAX_SIZE; i++) {
		if (data[i] == '\0') {
			return i;
		}
	}
	return StringStatsData::MAX_STRING_MINMAX_SIZE;
}

string StringStats::ToString(const BaseStatistics &stats) {
	auto &string_data = StringStats::GetDataUnsafe(stats);
	uint32_t min_len = GetValidMinMaxSubstring(string_data.min);
	uint32_t max_len = GetValidMinMaxSubstring(string_data.max);
	return StringUtil::Format("[Min: %s, Max: %s, Has Unicode: %s, Max String Length: %s]",
	                          Blob::ToString(string_t(const_char_ptr_cast(string_data.min), min_len)),
	                          Blob::ToString(string_t(const_char_ptr_cast(string_data.max), max_len)),
	                          string_data.has_unicode ? "true" : "false",
	                          string_data.has_max_string_length ? to_string(string_data.max_string_length) : "?");
}

void StringStats::Verify(const BaseStatistics &stats, Vector &vector, const SelectionVector &sel, idx_t count) {
	auto &string_data = StringStats::GetDataUnsafe(stats);

	UnifiedVectorFormat vdata;
	vector.ToUnifiedFormat(count, vdata);
	auto data = UnifiedVectorFormat::GetData<string_t>(vdata);
	for (idx_t i = 0; i < count; i++) {
		auto idx = sel.get_index(i);
		auto index = vdata.sel->get_index(idx);
		if (!vdata.validity.RowIsValid(index)) {
			continue;
		}
		auto value = data[index];
		auto data = value.GetData();
		auto len = value.GetSize();
		// LCOV_EXCL_START
		// 1. max string length check
		if (string_data.has_max_string_length && len > string_data.max_string_length) {
			throw InternalException(
			    "Statistics mismatch: string value exceeds maximum string length.\nStatistics: %s\nVector: %s",
			    stats.ToString(), vector.ToString(count));
		}
		// 2. unicode check
		if (stats.GetType().id() == LogicalTypeId::VARCHAR && !string_data.has_unicode) {
			auto unicode = Utf8Proc::Analyze(data, len);
			if (unicode == UnicodeType::UTF8) {
				throw InternalException("Statistics mismatch: string value contains unicode, but statistics says it "
				                        "shouldn't.\nStatistics: %s\nVector: %s",
				                        stats.ToString(), vector.ToString(count));
			} else if (unicode == UnicodeType::INVALID) {
				throw InternalException("Invalid unicode detected in vector: %s", vector.ToString(count));
			}
		}
		// 3. min/max zonemap check
		if (StringValueComparison(const_data_ptr_cast(data),
		                          MinValue<idx_t>(len, StringStatsData::MAX_STRING_MINMAX_SIZE), string_data.min) < 0) {
			throw InternalException("Statistics mismatch: value is smaller than min.\nStatistics: %s\nVector: %s",
			                        stats.ToString(), vector.ToString(count));
		}
		if (StringValueComparison(const_data_ptr_cast(data),
		                          MinValue<idx_t>(len, StringStatsData::MAX_STRING_MINMAX_SIZE), string_data.max) > 0) {
			throw InternalException("Statistics mismatch: value is bigger than max.\nStatistics: %s\nVector: %s",
			                        stats.ToString(), vector.ToString(count));
		}

		// 4) Prefix Bloom Filter check
		if (string_data.has_pbf) {
			for (int p = 0; p < StringStatsData::NUM_PREFIXES; p++) {
				uint32_t prefix_len = string_data.prefixes[p].level;
				if (prefix_len == 0 || len < prefix_len) {
					continue;
				}

				hash_t hash_val = Hash(reinterpret_cast<const char *>(data), prefix_len);
				uint32_t bit_index = hash_val % StringStatsData::NUM_BITS;
				uint32_t byte_idx = bit_index / 8;
				uint32_t bit_offset = bit_index % 8;

				uint8_t byte = string_data.prefixes[p].bits[byte_idx];
				if (((byte >> bit_offset) & 1) == 0) {
					throw InternalException(
					    "Statistics mismatch: prefix bloom filter bit not set for value.\nStatistics: %s\nVector: %s",
					    stats.ToString(), vector.ToString(count));
				}
			}
		}
		// LCOV_EXCL_STOP
	}
}

void StringStats::Init_PBF(StringStatsData& string_data){
	string_data.has_pbf = true;
	int level = 1;
	for(int i = 0; i < StringStatsData::NUM_PREFIXES; ++i){
		string_data.prefixes[i].level = level;
		for(int j = 0; j < StringStatsData::NUM_BYTES; ++j){
			string_data.prefixes[i].bits[j] = 0;
		}
		level *= 2;
	}
}

// Returns a PrefixQuery whose `prefixes` field contains the 1/2/4/8-byte
// prefixes derived from `constant` for comparison operators (=, >, <, >=, <=).
// Other expression types return an empty PrefixQuery.
// NOTE: BETWEEN is rewritten into two comparisons, so this function is invoked once per boundary and we compute prefixes for each side.
PrefixQuery StringStats::GetPrefixCandidates(ExpressionType comp_type,
                                             const std::string &constant) {
    PrefixQuery res;
    std::vector<std::string> prefixes;
    const std::vector<size_t> prefix_lengths = {1, 2, 4, 8};

    auto add_prefixes = [&]() {
        for (auto len : prefix_lengths) {
            if (len > constant.size()) break;
            prefixes.push_back(constant.substr(0, len));
        }
    };

    switch (comp_type) {
	// For filters like s = 'abc' or s IS NOT DISTINCT FROM 'abc'
    case ExpressionType::COMPARE_EQUAL:
    case ExpressionType::COMPARE_NOT_DISTINCT_FROM:
		std::cout << "[PBF]   case: EQUAL / NOT_DISTINCT_FROM\n";
        add_prefixes();
        break;

    default:
        break;
    }

    res.prefixes = std::move(prefixes);	
	return res;
}

} // namespace duckdb
