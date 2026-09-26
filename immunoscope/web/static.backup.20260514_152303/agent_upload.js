/**
 * ImmunoScope Agent with Trajectory Upload Support
 */

class TrajectoryManager {
    constructor() {
        this.trajectories = new Map();
        this.activeTrajectoryId = null;
        this.uploadEndpoint = '/api/agent/upload';

        this.setupUploadHandlers();
    }

    setupUploadHandlers() {
        const uploadArea = document.getElementById('uploadArea');
        const fileInput = document.getElementById('fileInput');

        // Click to upload
        uploadArea.addEventListener('click', () => fileInput.click());

        // File input change
        fileInput.addEventListener('change', (e) => {
            this.handleFiles(e.target.files);
            fileInput.value = ''; // Reset input
        });

        // Drag and drop
        uploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            uploadArea.classList.add('dragover');
        });

        uploadArea.addEventListener('dragleave', () => {
            uploadArea.classList.remove('dragover');
        });

        uploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('dragover');
            this.handleFiles(e.dataTransfer.files);
        });
    }

    async handleFiles(files) {
        const fileArray = Array.from(files);

        // Group files by trajectory (same base name)
        const groups = this.groupFilesByTrajectory(fileArray);

        for (const [baseName, fileList] of Object.entries(groups)) {
            const trajectoryId = this.generateId();
            const trajectory = {
                id: trajectoryId,
                name: baseName,
                files: fileList.map(f => ({
                    name: f.name,
                    size: f.size,
                    type: this.getFileType(f.name),
                    file: f
                })),
                status: 'uploading',
                uploadProgress: 0,
                serverPath: null
            };

            this.trajectories.set(trajectoryId, trajectory);
            this.renderTrajectory(trajectory);

            // Upload files
            await this.uploadTrajectory(trajectory);
        }
    }

    groupFilesByTrajectory(files) {
        const groups = {};

        for (const file of files) {
            // Extract base name (remove extension)
            const baseName = file.name.replace(/\.(xtc|tpr|pdb|gro|trr|dcd)$/i, '');

            if (!groups[baseName]) {
                groups[baseName] = [];
            }
            groups[baseName].push(file);
        }

        return groups;
    }

    async uploadTrajectory(trajectory) {
        const formData = new FormData();
        formData.append('trajectory_id', trajectory.id);
        formData.append('trajectory_name', trajectory.name);

        for (const fileInfo of trajectory.files) {
            formData.append('files', fileInfo.file);
        }

        try {
            const response = await fetch(this.uploadEndpoint, {
                method: 'POST',
                body: formData
            });

            if (!response.ok) {
                throw new Error(`Upload failed: ${response.statusText}`);
            }

            const result = await response.json();

            trajectory.status = 'ready';
            trajectory.serverPath = result.path;
            trajectory.uploadProgress = 100;

            this.updateTrajectoryUI(trajectory);

            // Auto-select first trajectory
            if (this.trajectories.size === 1) {
                this.setActiveTrajectory(trajectory.id);
            }

        } catch (error) {
            console.error('Upload error:', error);
            trajectory.status = 'error';
            trajectory.error = error.message;
            this.updateTrajectoryUI(trajectory);
        }
    }

    renderTrajectory(trajectory) {
        const listContainer = document.getElementById('trajectoryList');

        const item = document.createElement('div');
        item.className = 'trajectory-item';
        item.id = `trajectory-${trajectory.id}`;
        item.innerHTML = this.getTrajectoryHTML(trajectory);

        listContainer.appendChild(item);
    }

    updateTrajectoryUI(trajectory) {
        const item = document.getElementById(`trajectory-${trajectory.id}`);
        if (item) {
            item.innerHTML = this.getTrajectoryHTML(trajectory);

            // Re-attach event listeners
            this.attachTrajectoryListeners(trajectory.id);
        }
    }

    getTrajectoryHTML(trajectory) {
        const statusClass = trajectory.status;
        const statusText = {
            'uploading': '⏳ Uploading',
            'ready': '✓ Ready',
            'error': '✗ Error'
        }[trajectory.status];

        const filesHTML = trajectory.files.map(f =>
            `<div>📄 ${f.name} (${this.formatSize(f.size)})</div>`
        ).join('');

        const progressHTML = trajectory.status === 'uploading'
            ? `<div class="upload-progress">
                   <div class="upload-progress-bar" style="width: ${trajectory.uploadProgress}%"></div>
               </div>`
            : '';

        const actionsHTML = trajectory.status === 'ready'
            ? `<div class="trajectory-actions">
                   <button class="btn-small btn-analyze" onclick="trajectoryManager.analyzeTrajectory('${trajectory.id}')">
                       Analyze
                   </button>
                   <button class="btn-small btn-remove" onclick="trajectoryManager.removeTrajectory('${trajectory.id}')">
                       Remove
                   </button>
               </div>`
            : trajectory.status === 'error'
            ? `<div class="trajectory-actions">
                   <button class="btn-small btn-remove" onclick="trajectoryManager.removeTrajectory('${trajectory.id}')">
                       Remove
                   </button>
               </div>`
            : '';

        return `
            <div class="trajectory-header">
                <div class="trajectory-name">
                    🧬 ${trajectory.name}
                </div>
                <span class="trajectory-status ${statusClass}">${statusText}</span>
            </div>
            <div class="trajectory-files">
                ${filesHTML}
            </div>
            ${progressHTML}
            ${actionsHTML}
        `;
    }

    attachTrajectoryListeners(trajectoryId) {
        const item = document.getElementById(`trajectory-${trajectoryId}`);
        if (item) {
            item.onclick = (e) => {
                if (!e.target.closest('button')) {
                    this.setActiveTrajectory(trajectoryId);
                }
            };
        }
    }

    setActiveTrajectory(trajectoryId) {
        // Remove active class from all
        document.querySelectorAll('.trajectory-item').forEach(item => {
            item.classList.remove('active');
        });

        // Add active class to selected
        const item = document.getElementById(`trajectory-${trajectoryId}`);
        if (item) {
            item.classList.add('active');
            this.activeTrajectoryId = trajectoryId;
        }
    }

    analyzeTrajectory(trajectoryId) {
        const trajectory = this.trajectories.get(trajectoryId);
        if (!trajectory || trajectory.status !== 'ready') {
            return;
        }

        this.setActiveTrajectory(trajectoryId);

        // Auto-fill message input
        const input = document.getElementById('messageInput');
        input.value = `Analyze the trajectory "${trajectory.name}" at ${trajectory.serverPath}`;
        input.focus();
    }

    removeTrajectory(trajectoryId) {
        if (confirm('Remove this trajectory?')) {
            this.trajectories.delete(trajectoryId);
            const item = document.getElementById(`trajectory-${trajectoryId}`);
            if (item) {
                item.remove();
            }

            if (this.activeTrajectoryId === trajectoryId) {
                this.activeTrajectoryId = null;
            }
        }
    }

    getActiveTrajectory() {
        return this.activeTrajectoryId
            ? this.trajectories.get(this.activeTrajectoryId)
            : null;
    }

    getFileType(filename) {
        const ext = filename.split('.').pop().toLowerCase();
        const types = {
            'xtc': 'trajectory',
            'trr': 'trajectory',
            'dcd': 'trajectory',
            'tpr': 'topology',
            'pdb': 'structure',
            'gro': 'structure'
        };
        return types[ext] || 'unknown';
    }

    formatSize(bytes) {
        if (bytes < 1024) return bytes + ' B';
        if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
        if (bytes < 1024 * 1024 * 1024) return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
        return (bytes / (1024 * 1024 * 1024)).toFixed(1) + ' GB';
    }

    generateId() {
        return 'traj_' + Date.now() + '_' + Math.random().toString(36).substr(2, 9);
    }
}


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

        // Add trajectory context if available
        const activeTrajectory = trajectoryManager.getActiveTrajectory();
        let messageContent = content;

        if (activeTrajectory && activeTrajectory.status === 'ready') {
            // Include trajectory context in the message
            messageContent = content + `\n\n[Context: Using trajectory "${activeTrajectory.name}" at ${activeTrajectory.serverPath}]`;
        }

        this.addUserMessage(content, activeTrajectory);

        this.ws.send(JSON.stringify({
            type: 'message',
            content: messageContent,
            trajectory_context: activeTrajectory ? {
                id: activeTrajectory.id,
                name: activeTrajectory.name,
                path: activeTrajectory.serverPath,
                files: activeTrajectory.files.map(f => ({
                    name: f.name,
                    type: f.type
                }))
            } : null
        }));

        this.showTypingIndicator();
        this.disableInput();

        document.getElementById('messageInput').value = '';
        this.adjustTextareaHeight();
    }

    addUserMessage(content, trajectory) {
        const chatContainer = document.getElementById('chatContainer');

        const welcome = chatContainer.querySelector('.welcome-message');
        if (welcome) {
            welcome.remove();
        }

        const messageDiv = document.createElement('div');
        messageDiv.className = 'message user';

        const contextBadge = trajectory && trajectory.status === 'ready'
            ? `<div class="context-badge">📂 ${trajectory.name}</div>`
            : '';

        messageDiv.innerHTML = `
            <div class="message-content">
                ${this.escapeHtml(content)}
                ${contextBadge}
            </div>
            <div class="message-avatar">U</div>
        `;

        chatContainer.appendChild(messageDiv);
        this.scrollToBottom();
    }

    startAgentMessage() {
        if (this.currentMessageDiv) {
            return;
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

        input.addEventListener('input', () => this.adjustTextareaHeight());

        input.addEventListener('keydown', (e) => {
            if (e.key === 'Enter' && !e.shiftKey) {
                e.preventDefault();
                this.sendMessage(input.value);
            }
        });

        button.addEventListener('click', () => {
            this.sendMessage(input.value);
        });

        setInterval(() => {
            if (this.isConnected) {
                this.ws.send(JSON.stringify({ type: 'ping' }));
            }
        }, 30000);
    }
}

// Global instances
let trajectoryManager;
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
    trajectoryManager = new TrajectoryManager();
    agentClient = new AgentClient();
});
