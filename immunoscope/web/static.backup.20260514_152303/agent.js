/**
 * ImmunoScope Agent WebSocket Client
 */

class AgentClient {
    constructor() {
        this.ws = null;
        this.sessionId = null;
        this.isConnected = false;
        this.reconnectAttempts = 0;
        this.maxReconnectAttempts = 5;
        this.currentMessageDiv = null;

        this.init();
    }

    init() {
        this.connectWebSocket();
        this.setupEventListeners();
    }

    connectWebSocket() {
        const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
        const wsUrl = `${protocol}//${window.location.host}/api/agent/ws`;

        console.log('Connecting to:', wsUrl);
        this.ws = new WebSocket(wsUrl);

        this.ws.onopen = () => {
            console.log('WebSocket connected');
            this.isConnected = true;
            this.reconnectAttempts = 0;
            this.updateStatus('Connected', true);
            this.hideError();
        };

        this.ws.onmessage = (event) => {
            try {
                const message = JSON.parse(event.data);
                this.handleMessage(message);
            } catch (e) {
                console.error('Failed to parse message:', e);
            }
        };

        this.ws.onerror = (error) => {
            console.error('WebSocket error:', error);
            this.showError('Connection error. Please check your network.');
        };

        this.ws.onclose = () => {
            console.log('WebSocket closed');
            this.isConnected = false;
            this.updateStatus('Disconnected', false);

            // Attempt reconnection
            if (this.reconnectAttempts < this.maxReconnectAttempts) {
                this.reconnectAttempts++;
                const delay = Math.min(1000 * Math.pow(2, this.reconnectAttempts), 10000);
                console.log(`Reconnecting in ${delay}ms (attempt ${this.reconnectAttempts})`);
                setTimeout(() => this.connectWebSocket(), delay);
            } else {
                this.showError('Connection lost. Please refresh the page.');
            }
        };
    }

    handleMessage(message) {
        console.log('Received:', message);

        switch (message.type) {
            case 'session_info':
                this.sessionId = message.session_id;
                console.log('Session ID:', this.sessionId);
                break;

            case 'text_delta':
                this.appendAgentText(message.text);
                break;

            case 'tool_use_start':
                this.showToolCall(message.name, message.id, 'running');
                break;

            case 'tool_use_end':
                this.updateToolCall(message.id, message.ok ? 'success' : 'error', message.name);
                break;

            case 'turn_complete':
                this.finishAgentMessage();
                this.hideTypingIndicator();
                this.enableInput();
                break;

            case 'error':
                this.showError(message.message);
                this.hideTypingIndicator();
                this.enableInput();
                break;

            case 'aborted':
                this.appendAgentText('\n\n[Operation aborted by user]');
                this.finishAgentMessage();
                this.hideTypingIndicator();
                this.enableInput();
                break;

            case 'pong':
                // Heartbeat response
                break;

            default:
                console.warn('Unknown message type:', message.type);
        }
    }

    sendMessage(content) {
        if (!this.isConnected) {
            this.showError('Not connected to server');
            return;
        }

        if (!content.trim()) {
            return;
        }

        // Add user message to chat
        this.addUserMessage(content);

        // Send to server
        this.ws.send(JSON.stringify({
            type: 'message',
            content: content
        }));

        // Show typing indicator
        this.showTypingIndicator();
        this.disableInput();

        // Clear input
        document.getElementById('messageInput').value = '';
        this.adjustTextareaHeight();
    }

    addUserMessage(content) {
        const chatContainer = document.getElementById('chatContainer');

        // Remove welcome message if present
        const welcome = chatContainer.querySelector('.welcome-message');
        if (welcome) {
            welcome.remove();
        }

        const messageDiv = document.createElement('div');
        messageDiv.className = 'message user';
        messageDiv.innerHTML = `
            <div class="message-content">${this.escapeHtml(content)}</div>
            <div class="message-avatar">U</div>
        `;

        chatContainer.appendChild(messageDiv);
        this.scrollToBottom();
    }

    startAgentMessage() {
        if (this.currentMessageDiv) {
            return; // Already have an active message
        }

        const chatContainer = document.getElementById('chatContainer');

        this.currentMessageDiv = document.createElement('div');
        this.currentMessageDiv.className = 'message agent';
        this.currentMessageDiv.innerHTML = `
            <div class="message-avatar">A</div>
            <div class="message-content"></div>
        `;

        chatContainer.appendChild(this.currentMessageDiv);
        this.scrollToBottom();
    }

    appendAgentText(text) {
        if (!this.currentMessageDiv) {
            this.startAgentMessage();
        }

        const contentDiv = this.currentMessageDiv.querySelector('.message-content');
        contentDiv.textContent += text;
        this.scrollToBottom();
    }

    finishAgentMessage() {
        this.currentMessageDiv = null;
    }

    showToolCall(toolName, toolId, status = 'running') {
        if (!this.currentMessageDiv) {
            this.startAgentMessage();
        }

        const contentDiv = this.currentMessageDiv.querySelector('.message-content');

        const toolDiv = document.createElement('div');
        toolDiv.className = `tool-call ${status}`;
        toolDiv.id = `tool-${toolId}`;
        toolDiv.innerHTML = `
            <span class="tool-icon">${status === 'running' ? '⏳' : '🔧'}</span>
            <span>${toolName}</span>
        `;

        contentDiv.appendChild(toolDiv);
        this.scrollToBottom();
    }

    updateToolCall(toolId, status, toolName) {
        const toolDiv = document.getElementById(`tool-${toolId}`);
        if (toolDiv) {
            toolDiv.className = `tool-call ${status}`;
            const icon = status === 'success' ? '✓' : '✗';
            toolDiv.innerHTML = `
                <span class="tool-icon">${icon}</span>
                <span>${toolName}</span>
            `;
        }
    }

    showTypingIndicator() {
        const chatContainer = document.getElementById('chatContainer');

        let indicator = document.getElementById('typingIndicator');
        if (!indicator) {
            indicator = document.createElement('div');
            indicator.id = 'typingIndicator';
            indicator.className = 'message agent';
            indicator.innerHTML = `
                <div class="message-avatar">A</div>
                <div class="typing-indicator active">
                    <div class="typing-dots">
                        <span></span>
                        <span></span>
                        <span></span>
                    </div>
                </div>
            `;
            chatContainer.appendChild(indicator);
        }

        indicator.style.display = 'flex';
        this.scrollToBottom();
    }

    hideTypingIndicator() {
        const indicator = document.getElementById('typingIndicator');
        if (indicator) {
            indicator.remove();
        }
    }

    updateStatus(text, connected) {
        const statusText = document.getElementById('statusText');
        const statusDot = document.getElementById('statusDot');

        statusText.textContent = text;

        if (connected) {
            statusDot.classList.remove('disconnected');
        } else {
            statusDot.classList.add('disconnected');
        }
    }

    showError(message) {
        const errorBanner = document.getElementById('errorBanner');
        errorBanner.textContent = message;
        errorBanner.classList.add('active');

        // Auto-hide after 5 seconds
        setTimeout(() => this.hideError(), 5000);
    }

    hideError() {
        const errorBanner = document.getElementById('errorBanner');
        errorBanner.classList.remove('active');
    }

    enableInput() {
        const input = document.getElementById('messageInput');
        const button = document.getElementById('sendButton');
        input.disabled = false;
        button.disabled = false;
    }

    disableInput() {
        const input = document.getElementById('messageInput');
        const button = document.getElementById('sendButton');
        input.disabled = true;
        button.disabled = true;
    }

    scrollToBottom() {
        const chatContainer = document.getElementById('chatContainer');
        chatContainer.scrollTop = chatContainer.scrollHeight;
    }

    adjustTextareaHeight() {
        const textarea = document.getElementById('messageInput');
        textarea.style.height = 'auto';
        textarea.style.height = Math.min(textarea.scrollHeight, 150) + 'px';
    }

    escapeHtml(text) {
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }

    setupEventListeners() {
        const input = document.getElementById('messageInput');
        const button = document.getElementById('sendButton');

        // Auto-resize textarea
        input.addEventListener('input', () => this.adjustTextareaHeight());

        // Send on Enter (Shift+Enter for new line)
        input.addEventListener('keydown', (e) => {
            if (e.key === 'Enter' && !e.shiftKey) {
                e.preventDefault();
                this.sendMessage(input.value);
            }
        });

        // Send button click
        button.addEventListener('click', () => {
            this.sendMessage(input.value);
        });

        // Heartbeat
        setInterval(() => {
            if (this.isConnected) {
                this.ws.send(JSON.stringify({ type: 'ping' }));
            }
        }, 30000); // Every 30 seconds
    }
}

// Global functions
let agentClient;

function sendMessage() {
    const input = document.getElementById('messageInput');
    agentClient.sendMessage(input.value);
}

function sendExample(text) {
    const input = document.getElementById('messageInput');
    input.value = text;
    agentClient.sendMessage(text);
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    agentClient = new AgentClient();
});
