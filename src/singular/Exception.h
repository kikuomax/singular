#ifndef _SINGULAR_EXCEPTION_H
#define _SINGULAR_EXCEPTION_H

#include <string>

namespace singular {

	/** Exception. */
	class Exception {
	private:
		/** Message of this exception. */
		std::string message;
	public:
		/**
		 * Constructs with a message.
		 *
		 * @param message
		 *     Message of the exception.
		 */
		inline Exception (const std::string& message) : message(message) {}
	};

}

#endif
